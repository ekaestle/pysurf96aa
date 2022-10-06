"""
This script and the underlying Fortran scripts are based on
Computer Programs In Seismology (Hermann, 2013; https://www.eas.slu.edu/eqc/eqccps.html),
pysurf96 (https://opensourcelibs.com/lib/pysurf96) and
DAzimSurfTomo (Liu et al. 2019)

Liu, C., Yao, H., Yang, H.Y., Shen, W., Fang, H., Hu, S. and Qiao, L., 2019. Direct inversion for three‐dimensional shear wave speed azimuthal anisotropy based on surface wave ray tracing: Methodology and application to Yunnan, southwest China. Journal of Geophysical Research: Solid Earth, 124(11), pp.11394-11413.
"""

import numpy as num
from .surfdisp96aa_ext import surfdisp96, depthkernel  # noqa


MAXLAYER = 100
MAXPERIODS = 60


class Surf96Error(Exception):
    pass


def surf96aa(thickness, vp, vs, rho, anisoamp, anisodir, periods, nrefine=1,
             wave='rayleigh', mode=1,velocity='phase',flat_earth=False,return_sensitivities=False):

    if wave != 'rayleigh':
        raise Exception("anisotropic calculations currently only implemented for Rayleigh waves")
    if mode != 1:
        raise Exception("currently only fundamental mode calculations are supported (mode=1)")

    assert thickness.size == vp.size == vs.size == rho.size == anisoamp.size == anisodir.size, \
        'thickness, vp/vs velocities and rho have different sizes.'
    assert thickness.ndim == vp.ndim == vs.ndim == rho.ndim == 1, \
        'thickness, vp/vs velocities or rho have more than one dimension'
    assert thickness.size <= MAXLAYER, 'maximum number of layers is 100'
    assert nrefine*(thickness.size-1)+1 <= MAXLAYER, 'the number of layers times nrefine must not exceed 100'
    assert periods.size <= MAXPERIODS, 'maximum number of periods is 60'
    assert wave in ['love', 'rayleigh'], 'wave has to be love or rayleigh'
    assert velocity in ['group', 'phase'], 'velocity has to be group or phase'
    assert mode > 0, 'mode has to be at least 1 (fundamental mode)'

    iflsph = 1 if flat_earth else 0
    iwave = 1 if wave == 'love' else 2
    igr = 0 if velocity == 'phase' else 1
    mode = int(mode)

    nperiods = len(periods)
    nlayers = len(vp)

    _thk = num.empty(MAXLAYER)
    _vp = num.empty(MAXLAYER)
    _vs = num.empty(MAXLAYER)
    _rho = num.empty(MAXLAYER)

    _thk[:nlayers] = thickness
    _vp[:nlayers] = vp
    _vs[:nlayers] = vs
    _rho[:nlayers] = rho

    iflsph = 1 if flat_earth else 0
    iwave = 1 if wave == 'love' else 2
    igr = 0 if velocity == 'phase' else 1
    mode = int(mode)

    t = num.empty(MAXPERIODS)
    t[:nperiods] = periods

    c_iso = num.zeros(MAXPERIODS)
    
    # Lsen_Gsc is the first term in eq. (13) in Liu et al. 2019
    Lsen_Gsc,dcR_dA,dcR_dL,error = depthkernel(_thk, _vp, _vs, _rho,nlayers,nrefine,iflsph,iwave,
                        mode,igr,nperiods,t,c_iso)

    if error > 0:
        raise Surf96Error(
            'surf96 threw an error! '
            'This may be due to low velocity zone causing'
            ' reverse phase velocity dispersion,'
            ' and mode jumping. Due to looking for Love waves in a halfspace'
            ' which is OK if there are Rayleigh data.')


    nlayers_refined = Lsen_Gsc.shape[1]
    dvrel =   num.reshape(num.tile(num.repeat(anisoamp[:-1],nrefine),nperiods),(nperiods,nlayers_refined))

    # calculate dc according to eq. 13. (is sqrt((Gc/L)**2+(Gs/L)**2) = dvrel??)
    #fastdir = num.reshape(num.tile(num.repeat(anisodir[:-1],nrefine),nperiods),(nperiods,nlayers_refined))
    #angles=[0,20,40,60,80,90,110,130,150,170] # calculate dc for these propagation directions
    #if not hasattr(angles,'__iter__'):
    #    angles = [angles]
    #dc = num.zeros((nperiods,len(angles)))
    #for i,propdir in enumerate(angles):
    #    dc[:,i] = num.sum(Lsen_Gsc * dvrel*num.cos(2*(fastdir-propdir)),axis=1)

    # To calculate the velocity anisotropy we also need the terms Gc/L and Gs/L
    # Gc,s are the azimuthal variations of L, where L=rho*Vs
    # rho cancels out, dVs remains (is this correct?)
    Gc = dvrel * num.reshape(num.tile(num.repeat(num.cos(2*anisodir[:-1]),nrefine),nperiods),(nperiods,nlayers_refined))
    Gs = dvrel * num.reshape(num.tile(num.repeat(num.sin(2*anisodir[:-1]),nrefine),nperiods),(nperiods,nlayers_refined))
    # integration of eq. 13:
    C1 = num.sum(Lsen_Gsc * Gc, axis=1)
    C2 = num.sum(Lsen_Gsc * Gs, axis=1)
    # get anisotropic amplitude and angle (are the 0.5 factors correct? Seem to be necessary)
    c_aa_amp = 0.5*num.sqrt(C1**2 + C2**2)
    c_aa_ang = 0.5*num.arctan2(C2,C1)

    c_iso = c_iso[:nperiods]

    if return_sensitivities:
        return c_iso,c_aa_amp,c_aa_ang,Lsen_Gsc, dcR_dA, dcR_dL
    else:
        return c_iso,c_aa_amp,c_aa_ang


def surf96(thickness, vp, vs, rho, periods,
           wave='love', mode=1, velocity='group', flat_earth=True):
    '''Calculate synthetic surface wave dispersion curves

    A slim Fortran wrapper around surf96 from Computer Programs in Seismology
    from R. Hermann (2013).

    Parameters
    ----------
    thickness : numpy.array
        Layer thickness in kilometers
    vp : numpy.array
        Layer Vp velocity
    vs : numpy.array
        Layer Vs velocity
    rho : numpy.array
        Layer density in g/m^3
    periods : numpy.array
        The periods in seconds, where wave velocity is calculated
    wave : str
        The wave type, "love" or "rayleigh"
    mode : int
        Mode of the wave, 1: fundamental, 2: second-mode, etc...
    velocity : str
        "group" or "phase" velocity
    flat_earth : bool
        Assume a flat earth

    Returns
    -------
    numpy.array
        The surface wave velocities at defined periods.

    Raises
    ------
    Surf96Error
        If surf96 fortran code raises an error,
        this may be due to low velocity zone.
    '''
    assert thickness.size == vp.size == vs.size == rho.size, \
        'thickness, vp/vs velocities and rho have different sizes.'
    assert thickness.ndim == vp.ndim == vs.ndim == rho.ndim == 1, \
        'thickness, vp/vs velocities or rho have more than one dimension'
    assert thickness.size <= MAXLAYER, 'maximum number of layers is 100'
    assert periods.size <= MAXPERIODS, 'maximum number of periods is 60'
    assert wave in ['love', 'rayleigh'], 'wave has to be love or rayleigh'
    assert velocity in ['group', 'phase'], 'velocity has to be group or phase'
    assert mode > 0, 'mode has to be at least 1 (fundamental mode)'

    nlayers = thickness.size
    kmax = periods.size

    _thk = num.empty(MAXLAYER)
    _vp = num.empty(MAXLAYER)
    _vs = num.empty(MAXLAYER)
    _rho = num.empty(MAXLAYER)

    _thk[:nlayers] = thickness
    _vp[:nlayers] = vp
    _vs[:nlayers] = vs
    _rho[:nlayers] = rho

    iflsph = 1 if flat_earth else 0
    iwave = 1 if wave == 'love' else 2
    igr = 0 if velocity == 'phase' else 1
    mode = int(mode)

    t = num.empty(MAXPERIODS)
    t[:kmax] = periods

    result = num.zeros(MAXPERIODS)

    error = surfdisp96(_thk, _vp, _vs, _rho, nlayers, iflsph, iwave,
                       mode, igr, kmax, t, result)
    if error > 0:
        raise Surf96Error(
            'surf96 threw an error! '
            'This may be due to low velocity zone causing'
            ' reverse phase velocity dispersion,'
            ' and mode jumping. Due to looking for Love waves in a halfspace'
            ' which is OK if there are Rayleigh data.')

    return result[:kmax]


def layermod2depthmod(thickness,parameters):

    if type(parameters)!=type(()):
        raise Exception("parameters has to be a tuple (vp,vs,...)")

    depth = num.append(0,num.repeat(num.cumsum(thickness),2)[:-1])
    depth[-1] = depth[-2]*1.1 # halfspace
    layerparams = num.repeat(parameters,2,axis=1)

    return num.vstack((depth,layerparams))


__all__ = ['surf96', 'surf96aa', 'layermod2depthmod','Surf96Error']

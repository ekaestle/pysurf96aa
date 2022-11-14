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
             wave='rayleigh', mode=1,velocity='phase',flat_earth=False,return_sensitivities=False,version=1):
    '''Calculate synthetic surface wave dispersion curves including 2Psi azimuthal anisotropy

    A Fortran wrapper around surf96 and tregn96 from Computer Programs in Seismology from R. Hermann (2013).

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
    anisoamp : numpy.array
        Layer anisotropic amplitude (a2) in relative units as in V = Viso * (1 + a2*cos(2*Psi2))
    anisodir : numpy.array
        Layer anisotropic directions (Psi2) in radians
    periods : numpy.array
        The periods in seconds, where wave velocity is calculated
    nrefine : integer
        subdivide layers into smaller layers, will only influence how
        the sensitivity kernels are sampled, otherwise identical results.
    wave : str
        The wave type, "love" or "rayleigh". Only Rayleigh waves are supported currently.
    mode : int
        Mode of the wave, 1: fundamental, 2: second-mode, etc... Only fundamental mode is supported currently.
    velocity : str
        "group" or "phase" velocity. Only phase velocities are supported currently.
    flat_earth : bool
        Assume a flat earth
    return_sensitivities : bool
        Returns the sensitivity kernels for the azimuthal anisotropy (dC/dA and dC/dL) as well as the joint,
        effective sensitivity as in Liu et al. 2019 eq. 13 (dC/dA*A+dC/dL*L).
        Returns tuple: c_iso,c_aa_amp,c_aa_ang,Lsen_Gsc, dcR_dA, dcR_dL
    version : integer
        Can be either 1 or 2. If 1, calculations are according to Liu et al. 2019 otherwise the method of
        Bodin et al. 2016 is applied. The results are very similar.

    Returns
    -------
    tuple of numpy.array
        c_iso: The surface wave velocities at defined periods.
        c_aa_amp: The relative azimuthal anisotropy amplitude as in c = c_iso * (1 + c_aa_amp * cos(2*c_aa_ang))
        c_aa_ang: The azimuthal anisotropy directions
        Lsen_Gsc (optional): The effective sensitivity kernels (size: periods*(layers-1))
        dcR_dA (optional): Partial derivative, i.e. the sensitivity to shallow P velocity (same size as Lsen_Gsc)
        dcR_dL (optional): Partial derivative, i.e. the sensitivity to deeper S velocity (same size as Lsen_Gsc)

    Raises
    ------
    Surf96Error
        If surf96 fortran code raises an error,
        this may be due to a low velocity zone.
    '''

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
    
    # dcR_dA and dcR_dL are the partial derivatives in Liu et al. 2019 (eq. 12) or Bodin et al. 2016 (eq A3,A4)
    Lsen_Gsc,dcR_dA,dcR_dL,error = depthkernel(_thk, _vp, _vs, _rho,nlayers,nrefine,iflsph,iwave,
                                               mode,igr,nperiods,t,c_iso)

    if error > 0:
        raise Surf96Error(
            'surf96 threw an error! '
            'This may be due to low velocity zone causing'
            ' reverse phase velocity dispersion,'
            ' and mode jumping. Due to looking for Love waves in a halfspace'
            ' which is OK if there are Rayleigh data.')

    nlayers_refined = dcR_dA.shape[1]
    # A = rho*vp**2
    # L = rho*vs**2
    # G = rho*vs*dVs (from Bodin et al. 2016, A5)
    if version==1: # as described in Liu et al. 2019 (eq. 11-13)
        # assumption: B/A = G/L
        # C1,2 = integral ( dcR_dA * B + dcR_dL * G ) = integral ( (G/L) * ( dcR_dA * A + dcR_dL * L ) ) = integral( G/L * Lsen_Gc )
        # no need to calculate Lsen_Gsc since it is calculated in the Fortran script directly
        #Lsen_Gsc = (dcR_dA*num.reshape(num.tile(num.repeat(rho[:-1]*vp[:-1]**2,nrefine),nperiods),(nperiods,nlayers_refined)) +
        #            dcR_dL*num.reshape(num.tile(num.repeat(rho[:-1]*vs[:-1]**2,nrefine),nperiods),(nperiods,nlayers_refined)) )
        # C1 = integral( G/L * Lsen_Gc ) = integral( G/L * Lsen_Gsc * cosine(2*anisodir) )
        # G/L = dVs/vs = 2*anisoamp (dVs is the peak to peak amplitude while anisoamp is the relative amplitude from V = Viso*(1+anisoamp*cos(2*anisodir)))
        C1 = num.sum(Lsen_Gsc*num.reshape(num.tile(num.repeat(2*anisoamp[:-1]*num.cos(2*anisodir[:-1]),nrefine),nperiods),(nperiods,nlayers_refined)),axis=1)
        # C2 = integral( G/L * Lsen_Gs ) = integral( G/L * Lsen_Gsc * sine(2*anisodir) )
        C2 = num.sum(Lsen_Gsc*num.reshape(num.tile(num.repeat(2*anisoamp[:-1]*num.sin(2*anisodir[:-1]),nrefine),nperiods),(nperiods,nlayers_refined)),axis=1)
    else: # as described in Bodin et al. 2016 (appendix)
        # assumption: Vp anisotropy points in the same direction and has an amplitude 1.5 times stronger than Vs (Obreski et al. 2010, Bodin et al. 2016)
        # if we set the scaling factor to 1, the results of both versions are the same.
        dVs = 2*anisoamp*vs # peak to peak absolute velocity deviation
        dVp = 2*1.5*anisoamp*vp # see assumption
        # C1 = integral ( dcR_dA * Bc + dcR_dL * Gc ) = C1 = integral ( dcR_dA * B * costerm + dcR_dL * G * costerm )
        # C2 = integral ( dcR_dA * Bs + dcR_dL * Gs ) = C1 = integral ( dcR_dA * B * sinterm + dcR_dL * G * sinterm )
        costerm = num.reshape(num.tile(num.repeat(num.cos(2*anisodir[:-1]),nrefine),nperiods),(nperiods,nlayers_refined))
        sinterm = num.reshape(num.tile(num.repeat(num.sin(2*anisodir[:-1]),nrefine),nperiods),(nperiods,nlayers_refined))
        G = num.reshape(num.tile(num.repeat(rho[:-1]*vs[:-1]*dVs[:-1],nrefine),nperiods),(nperiods,nlayers_refined))
        B = num.reshape(num.tile(num.repeat(rho[:-1]*vp[:-1]*dVp[:-1],nrefine),nperiods),(nperiods,nlayers_refined))
        # integration according to eq. A3, A4 of Bodin et al. 2016
        C1 = num.sum(dcR_dA*B*costerm + dcR_dL*G*costerm,axis=1)
        C2 = num.sum(dcR_dA*B*sinterm + dcR_dL*G*sinterm,axis=1)
    
    # get anisotropic amplitude and angle
    c_aa_amp = num.sqrt(C1**2 + C2**2) # absolute dc phase velocity deviation
    c_aa_ang = 0.5*num.arctan2(C2,C1)

    c_iso = c_iso[:nperiods]
    c_aa_amp /= c_iso # relative amplitude

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

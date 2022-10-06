# PySurf96AA
![Python23](https://img.shields.io/badge/python-2.7%20%7C%203.x-brightgreen.svg)

**This is based on PySurf96 (https://github.com/miili/pysurf96) and some of the Fortran routines are taken from Chuanming Liu (https://github.com/Chuanming-Liu/DAzimSurfTomo/).**

Next to the normal functionality of surf96, it includes the possibility to calculate azimuthally anisotropic Rayleigh phase dispersion. The current version may still contain errors, please check carefully before using the code.

_Modelling Surface Wave Dispersion Curves_

This is a slim wrapper around the program `surf96` from _Computer programs in seismology_ by R. Hermann (http://www.eas.slu.edu/eqc/eqccps.html) for forward modelling of Rayleigh and Love wave dispersion curves.

In this realisation the Fortran77 code is wrapped by `f2py`, which makes the forward computation approximately **8x faster** compared over calling a Python subprocess.

More useful software for seismology at https://pyrocko.org.

## Installation

This package is for Python 2 and Python 3.

Prerequisits are numpy and a Fortran77 compiler, like GNU GCC.

```
sudo python setup.py install
```

Or through pip:

```
pip install git+https://github.com/ekaestle/pysurf96aa
```

## Input parameters

```python
c_iso,aa_amp,aa_dir = surf96aa(thickness, vp, vs, rho, anisoamp, anisodir, 
    periods, nrefine=1, wave='rayleigh', mode=1, velocity='phase',
    flat_earth=False, return_sensitivities=False)
c_iso,aa_amp,aa_dir,Lsen,dcda,dcdl = surf96aa(thickness, vp, vs, rho, anisoamp, anisodir,
    periods, nrefine=1, wave='rayleigh', mode=1, velocity='phase',
    flat_earth=False, return_sensitivities=True)
```
**Parameter explanations**

thickness: thickness of each layer in km (last layer is halfspace, thickness is ignored)

vp: Vp velocity of each layer in km/s, same length as thickness

vs: Vs velocity of each layer in km/s

rho: density of each layer in g/cm³

anisoamp: relative anisotropic amplitude in each layer (0.05 = 5% anisotropy)

anisodir: fast axis direction in radians in each layer

periods: periods at which to return the phase dispersion

nrefine: split layers in _nrefine_ sublayers. Has no effect on the result, but gives smoother sensitivity kernels

wave: only 'rayleigh' is currently supported for azimuthal anisotropy calculations

mode: only fundamental mode (=1) is currently supported for azimuthal anisotropy calculations

velocity: only 'phase' is currently supported for azimuthal anisotropy calculations

flat_earth: apply flat earth transform

return_sensitivity: if _True_, the partial derivatives dC/dA and dC/dL as well as the total sensitivity as in the first term of eq. (13) of Liu et al. (2019) is returned.


## Example

**Isotropic**

This is the same as PySurf96
```python
import numpy as np
from pysurf96aa import surf96

# Define the velocity model in km and km/s
thickness = np.array([5., 23., 8., 0])
vs = np.array([2, 3.6, 3.8, 3.3])
vp = vs * 1.73
rho = vp * .32 + .77

# Periods we are interested in
periods = np.linspace(1., 20., 20)

velocities = surf96(thickness, vp, vs, rho, periods,
                    wave='rayleigh', mode=1, velocity='phase', flat_earth=False)
```

**Anisotropic**
```python
import numpy as np
from pysurf96aa import surf96aa

# Define the velocity model in km and km/s
thickness = np.array([5., 23., 8., 0])
vs = np.array([2, 3.6, 3.8, 3.3])
vp = vs * 1.73
rho = vp * .32 + .77
# azimuthal anisotropy amplitude (relative to vs)
aniso_amp = np.array([0.03,0.05,0.02,0.02])
# azimuthal anisotropy fast axis direction (in radians)
aniso_dir = np.array([30,-30,30.,0.])

# Periods we are interested in
periods = np.linspace(1., 20., 20)

c_iso,aa_amp,aa_dir = surf96aa(thickness, vp, vs, rho, aniso_amp, aniso_dir, periods,
                               wave='rayleigh', mode=1,velocity='phase',flat_earth=False)

```

There is also a _testrun.py_ script in the tutorial folder that can be used to test the functionalities and create some illustrative plots.


## Citations and Acknowledgments

> Herrmann, R. B. (2013) Computer programs in seismology: An evolving tool for instruction and research, Seism. Res. Lettr. 84, 1081-1088, doi:10.1785/0220110096

> Liu, C., Yao, H., Yang, H.Y., Shen, W., Fang, H., Hu, S. and Qiao, L., 2019. Direct inversion for three‐dimensional shear wave speed azimuthal anisotropy based on surface wave ray tracing: Methodology and application to Yunnan, southwest China. Journal of Geophysical Research: Solid Earth, 124(11), pp.11394-11413.

Thanks to Hongjian Fang for creating the Fortran subroutine (https://github.com/caiweicaiwei/SurfTomo)


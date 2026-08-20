"""Python sidecar of the Fortran sed_astrodust SED solver.

Exposes single-grain `sed_equilibrium` and `sed_stochastic` for
prototyping and verification. Use the Fortran library
(`sed/calc_sed.x astrodust`) when running the full size-distribution sum.

`sedust_h5` reads the optics products themselves -- the wavelength axis, the
cross-section tables and the extinction curve of `data/<model>/sedust_<model>.h5` --
and is the counterpart of `sed/src/sedust_product.f90`. It needs h5py; the
rest of this package does not, so it is imported lazily and a tree without
h5py keeps working.
"""
from .sed_from_cabs import (
    planck_lambda,
    equilibrium_T,
    sed_equilibrium,
    calc_P,
    sed_stochastic,
)
from .radiation_fields import J_Mathis

__all__ = [
    'planck_lambda',
    'equilibrium_T',
    'sed_equilibrium',
    'calc_P',
    'sed_stochastic',
    'J_Mathis',
    # sedust_h5
    'model_file',
    'read_grid',
    'read_kext',
    'read_qtable',
    'read_polarized',
    'describe',
]


def __getattr__(name):
    if name in ('model_file', 'read_grid', 'read_kext', 'read_qtable',
                'read_polarized', 'describe', 'components',
                'Grid', 'Kext', 'QTable', 'Polarized'):
        from . import sedust_h5
        return getattr(sedust_h5, name)
    raise AttributeError(f'module {__name__!r} has no attribute {name!r}')

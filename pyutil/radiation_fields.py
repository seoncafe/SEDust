"""Radiation-field constructors for the Python SED utility.

Mirrors sed/src/radfield.f90 :: J_Mathis (same piecewise
definition plus a CMB blackbody component) so callers can build a
Mathis-ISRF cell without re-deriving the formula. The default matches the
Fortran default ``use_mathis_corrected = .true.`` (Draine's corrected
4000 K dilution factor 1.65e-13 and 2.725 K CMB); passing
``corrected=False`` reproduces the literal Mathis 1983 values
(1.0e-13 / 2.9 K), which is what the Fortran CLI option 'mathis_orig'
selects.
"""
from __future__ import annotations
import numpy as np
try:
    from .sed_from_cabs import planck_lambda
except ImportError:                 # allow running as a plain script
    from sed_from_cabs import planck_lambda


def J_Mathis(U: float, lam_um: np.ndarray, corrected: bool = True) -> np.ndarray:
    """Mathis 1983 ISRF, scaled by intensity factor U, plus a CMB blackbody.

    Mirrors sed/src/radfield.f90 :: J_Mathis. The default
    ``corrected=True`` reproduces the Fortran default
    ``use_mathis_corrected = .true.``: Draine's 2008.02.02 correction of
    the 4000 K dilution factor (1.0e-13 -> 1.65e-13) together with the
    modern CMB temperature (2.9 -> 2.725 K). ``corrected=False`` recovers
    the literal Mathis 1983 values (1.0e-13 / 2.9 K), which is what the
    Fortran CLI option 'mathis_orig' selects. Only the 4000 K dilution
    factor and the CMB temperature differ between the two branches; the
    7500 K (1e-14) and 3000 K (4e-13) terms are identical in both. The CMB
    is added unscaled by U (matches the Fortran and the HD23 convention).

    Returns J_lambda in SI W m^-3 sr^-1 (the same convention as
    `radfield :: bbody` and `sed_astrodust_mod :: J_lam`).
    """
    if corrected:
        w_4000 = 1.65e-13
        T_cmb = 2.725
    else:
        w_4000 = 1.0e-13
        T_cmb = 2.9
    lam = np.asarray(lam_um, dtype=np.float64)
    J = np.zeros_like(lam)
    # Mathis piecewise (lambda in um)
    m1 = lam < 0.0912
    m2 = (lam >= 0.0912) & (lam < 0.110)
    m3 = (lam >= 0.110) & (lam < 0.134)
    m4 = (lam >= 0.134) & (lam < 0.250)
    m5 = lam >= 0.250
    J[m1] = 0.0
    J[m2] = 3069.0 * lam[m2]**3.4172
    J[m3] = 1.627
    J[m4] = 0.0566 * lam[m4]**(-1.6678)
    J[m5] = (1.0e-14 * planck_lambda(lam[m5], 7500.0)
             + w_4000 * planck_lambda(lam[m5], 4000.0)
             + 4.0e-13 * planck_lambda(lam[m5], 3000.0))
    # Scale by U then add the CMB
    J *= U
    J += planck_lambda(lam, T_cmb)
    return J

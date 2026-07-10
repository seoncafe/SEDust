"""Radiation-field constructors for the Python SED utility.

Mirrors astrodust/sed/src/radfield.f90 :: J_Mathis (and adds the same
piecewise definition + 2.9 K CMB component) so callers can build a
Mathis-ISRF cell without re-deriving the formula.
"""
from __future__ import annotations
import numpy as np
try:
    from .sed_from_cabs import planck_lambda
except ImportError:                 # allow running as a plain script
    from sed_from_cabs import planck_lambda


def J_Mathis(U: float, lam_um: np.ndarray) -> np.ndarray:
    """Mathis 1983 ISRF, scaled by intensity factor U, plus 2.9 K CMB.

    Returns J_lambda in SI W m^-3 sr^-1 (the same convention as
    `radfield :: bbody` and `sed_astrodust_mod :: J_lam`).
    """
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
             + 1.0e-13 * planck_lambda(lam[m5], 4000.0)
             + 4.0e-13 * planck_lambda(lam[m5], 3000.0))
    # Scale by U then add the CMB
    J *= U
    J += planck_lambda(lam, 2.9)
    return J

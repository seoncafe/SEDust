"""SED of a single dust grain from its absorption cross section and the
incident radiation field.

Two solve modes:

  equilibrium_T(lam_um, Cabs, J_lam) -> T_eq
  sed_equilibrium(lam_um, Cabs, J_lam) -> T_eq, lambda*I_lambda
      Solve absorbed-power = emitted-power for an equilibrium grain
      temperature; emission is C_abs * B_lambda(T_eq).

  sed_stochastic(lam_um, Cabs, J_lam, T_grid, U_grain) -> P_T, lambda*I_lambda
      Guhathakurta & Draine (1989) matrix solver for the temperature
      distribution P(T) of a stochastically heated grain, with the user-
      provided enthalpy U(T) [erg]. Emission is sum_i P(T_i) * C_abs *
      B_lambda(T_i).

This module mirrors the conventions of the Fortran sed_astrodust_mod
shipped with this repository: lambda in microns, C_abs in cm^2,
J_lambda in erg s^-1 cm^-2 um^-1 sr^-1 (Mathis-ISRF compatible), and
the returned lambda * I_lambda is in erg s^-1 sr^-1 per unit grain
(the caller multiplies by `dn_grain / dN_H` to get the emission per H
atom that HD23 reports).

Why a Python copy: this is the "downstream" form of the dust SED
solver — the function signature takes only what a 3D-RT cell needs to
provide (`lam`, `J_lam`, `C_abs`) and returns what the cell needs to
consume (`lam * I_lam`). Writing it in Python first makes it easy to
prototype radiation field shapes, plot intermediate quantities, and
validate the Fortran library on representative grains without the
full T-matrix + size-distribution pipeline.
"""

from __future__ import annotations
import numpy as np
from scipy.constants import h, c, k
from scipy.optimize import brentq
from scipy.integrate import simpson
from typing import Tuple


# ---------------------------------------------------------------------------
# Planck function
# ---------------------------------------------------------------------------

def planck_lambda(lam_um: np.ndarray, T: float) -> np.ndarray:
    """B_lambda(T) at lam_um, in W m^-3 sr^-1 (SI, per unit wavelength: W m^-2 m^-1 sr^-1).

    Mirrors the Fortran `radfield :: bbody` (SI). Multiply by 10 to get CGS
    (erg s^-1 cm^-3 sr^-1 per cm wavelength).
    """
    lam_m = np.asarray(lam_um, dtype=np.float64) * 1.0e-6
    x = h * c / (lam_m * k * T)
    # Use np.expm1 to keep precision in the long-wavelength (small x) limit.
    return 2.0 * h * c**2 / lam_m**5 / np.expm1(x)


# ---------------------------------------------------------------------------
# Equilibrium temperature
# ---------------------------------------------------------------------------

def absorbed_power(lam_um: np.ndarray, Cabs: np.ndarray, J_lam: np.ndarray) -> float:
    """Total absorbed power per grain (mixed-unit: cm^2 * SI Bλ * dλ).

    Multiply by 4*pi to get total flux integrated over solid angle.
    The kappJ value passed to the equilibrium balance is just the
    integrand of `int Cabs(lam) J_lam(lam) dlam`.
    """
    return simpson(Cabs * J_lam, x=lam_um)


def emitted_power(lam_um: np.ndarray, Cabs: np.ndarray, T: float) -> float:
    """`int Cabs(lam) B_lam(T) dlam` (mixed-unit: cm^2 * SI Bλ * dλ).

    Multiply by 4*pi for the total radiated luminosity of the grain.
    """
    return simpson(Cabs * planck_lambda(lam_um, T), x=lam_um)


def equilibrium_T(lam_um: np.ndarray,
                  Cabs: np.ndarray,
                  J_lam: np.ndarray,
                  T_lo: float = 2.7,
                  T_hi: float = 3.0e3) -> float:
    """Solve `int Cabs J_lam dlam = int Cabs B_lam(Teq) dlam` for Teq.

    The 4*pi cancels on both sides (it's the same on absorption and
    emission for an isotropic grain in LTE). brentq bracketing.
    """
    P_abs = absorbed_power(lam_um, Cabs, J_lam)
    if P_abs <= 0.0:
        return T_lo
    return brentq(lambda T: emitted_power(lam_um, Cabs, T) - P_abs, T_lo, T_hi)


def sed_equilibrium(lam_um: np.ndarray,
                    Cabs: np.ndarray,
                    J_lam: np.ndarray,
                    T_lo: float = 2.7,
                    T_hi: float = 3.0e3) -> Tuple[float, np.ndarray]:
    """Single-grain equilibrium SED.

    Returns
    -------
    T_eq : float
        Equilibrium temperature [K].
    lamI_lam : np.ndarray
        lambda * I_lambda in erg/s/sr per grain on the input lam_um grid.
        Convert to units per H atom by multiplying by `dn_grain/dN_H`.

    The 1e-3 factor converts the mixed-unit product
    `lam_um * Cabs[cm^2] * B_lam[SI W/m^3/sr]` into `erg/s/sr per grain`
    (see sed_astrodust_mod :: sed_solve, same formula).
    """
    Teq = equilibrium_T(lam_um, Cabs, J_lam, T_lo, T_hi)
    lamI_lam = lam_um * Cabs * planck_lambda(lam_um, Teq) * 1.0e-3
    return Teq, lamI_lam


# ---------------------------------------------------------------------------
# Stochastic heating (Guhathakurta & Draine 1989 matrix solver)
# ---------------------------------------------------------------------------

# Physical constants in CGS used by the matrix construction
_PI = np.pi
_FOURPI = 4.0 * _PI
_HC_ERG_UM = (h * c / 1e-6) * 1e7   # h*c in erg*um (lambda in um -> energy in erg)


def _build_kappB(lam_um: np.ndarray,
                 Cabs: np.ndarray,
                 T_grid: np.ndarray) -> np.ndarray:
    """`int Cabs(lam) B_lam(T) dlam` for each T in T_grid.

    Used by both calc_Teq and calc_P. Returns array of length len(T_grid)
    in mixed units (cm^2 * SI Bλ * dλ-um).
    """
    return np.array([emitted_power(lam_um, Cabs, T) for T in T_grid])


def calc_P(lam_um: np.ndarray,
           Cabs: np.ndarray,
           J_lam: np.ndarray,
           T_grid: np.ndarray,
           U_grain: np.ndarray,
           kappCMB: float = 0.0) -> np.ndarray:
    """Guhathakurta & Draine (1989) matrix solver for P(T).

    Mirrors `p_sub :: calc_P` from this repository (LINEAR_P branch).

    Parameters
    ----------
    lam_um   : (NW,) wavelength grid [um]
    Cabs     : (NW,) absorption cross section [cm^2] for the grain
    J_lam    : (NW,) mean intensity on the same lam_um grid [SI W/m^3/sr]
    T_grid   : (NT,) temperature grid [K], log-spaced and ideally already
               narrowed around the expected Teq for this grain
    U_grain  : (NT,) enthalpy U(T_grid) [erg]
    kappCMB  : scalar, 2.9 K Planck integral (subtracted from kappB to
               null out the heat-bath term in the absence of true CMB
               photons in the lam grid). Usually small; pass 0 to skip.

    Returns
    -------
    P : (NT,) probability density on T_grid (sums to 1).
    """
    NW = len(lam_um)
    NT = len(T_grid)
    assert Cabs.shape == (NW,)
    assert J_lam.shape == (NW,)
    assert U_grain.shape == (NT,)

    Amat = np.zeros((NT, NT))

    # Same kappB as in equilibrium calc, evaluated at every T grid point.
    kappB = _build_kappB(lam_um, Cabs, T_grid)

    # Down-transition rates (i+1 -> i): cooling by emitting kappB photons.
    H = U_grain
    for i1 in range(1, NT):
        kappB_tot = kappB[i1] - kappCMB
        Amat[i1-1, i1] = _FOURPI / (H[i1] - H[i1-1]) * kappB_tot

    # Up-transition rates (i1 -> i2 with i2 > i1): absorption of one
    # photon of energy H(i2) - H(i1).
    log_lam = np.log(lam_um)
    for i1 in range(NT-1):
        for i2 in range(i1+1, NT):
            ener = H[i2] - H[i1]
            if i2 == NT-1:
                delH = H[i2] * 0.5 * np.log(H[i2]/H[i2-1])
            else:
                delH = H[i2] * 0.5 * np.log(H[i2+1]/H[i2-1])
            wavl = _HC_ERG_UM / ener
            Jw = float(np.interp(np.log(wavl), log_lam, J_lam))
            Cx = float(np.interp(np.log(wavl), log_lam, Cabs))
            Amat[i2, i1] = _FOURPI * Cx * _HC_ERG_UM * delH / ener**3 * Jw

            # Highest bin: add transitions from photons more energetic
            # than what fits in this bin (GD89 closure).
            if i2 == NT-1:
                nwav = 51
                dlnwav = np.log(wavl/lam_um[0]) / (nwav - 1)
                wgrid = lam_um[0] * np.exp(np.arange(nwav) * dlnwav)
                Jwg = np.interp(np.log(wgrid), log_lam, J_lam)
                Cxg = np.interp(np.log(wgrid), log_lam, Cabs)
                Amat[i2, i1] += np.sum(_FOURPI * Cxg * wgrid**2 / _HC_ERG_UM * dlnwav * Jwg)

    # Cumulate: Bmat(j, i) = sum_{k>=j} Amat(k, i).
    # The solver is the same forward sweep as in p_sub.f90 LINEAR_P branch.
    Bmat = np.zeros((NT, NT))
    Bmat[NT-1, :NT-1] = Amat[NT-1, :NT-1]
    sumB = np.sum(Bmat[NT-1, :NT-1])
    for i1 in range(NT-1):
        for i2 in range(NT-2, i1, -1):
            Bmat[i2, i1] = Bmat[i2+1, i1] + Amat[i2, i1]
            sumB += Bmat[i2, i1]

    P = np.zeros(NT)
    if sumB > 0.0:
        P[0] = 1.0
        for j in range(1, NT):
            if Amat[j-1, j] > 0.0:
                P[j] = np.dot(Bmat[j, :j], P[:j]) / Amat[j-1, j]
            if P[j] > 1e50:
                P[:j+1] /= P[j]      # rescale to avoid overflow
        sP = P.sum()
        if sP > 0.0:
            P /= sP
    return P


def sed_stochastic(lam_um: np.ndarray,
                   Cabs: np.ndarray,
                   J_lam: np.ndarray,
                   T_grid: np.ndarray,
                   U_grain: np.ndarray,
                   kappCMB: float = 0.0) -> Tuple[np.ndarray, np.ndarray]:
    """Single-grain stochastic SED.

    Returns
    -------
    P : (NT,)
        Temperature probability density on T_grid (sums to 1).
    lamI_lam : (NW,)
        lambda * I_lambda in erg/s/sr per grain.
        sum over T:  sum_i P_i * lam * Cabs * B_lam(T_i) * 1e-3
    """
    P = calc_P(lam_um, Cabs, J_lam, T_grid, U_grain, kappCMB=kappCMB)
    lamI_lam = np.zeros_like(lam_um)
    for i, Ti in enumerate(T_grid):
        if P[i] > 0.0:
            lamI_lam += P[i] * Cabs * planck_lambda(lam_um, Ti)
    lamI_lam *= lam_um * 1.0e-3
    return P, lamI_lam

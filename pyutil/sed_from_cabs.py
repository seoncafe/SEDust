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
      provided enthalpy U(T) [erg] and the finite-bin upward transition
      rates of Draine & Li (2001). Emission is sum_i P(T_i) * C_abs *
      B_lambda(T_i).

Units follow the Fortran sed_astrodust_mod of this repository: lambda in
microns, C_abs in cm^2, and J_lambda in SI W m^-3 sr^-1 (Mathis-ISRF
compatible, the same convention as `radfield :: bbody` and `planck_lambda`
below). The returned lambda * I_lambda is in erg s^-1 sr^-1 per unit grain
(the caller multiplies by `dn_grain / dN_H` to get the emission per H atom
that HD23 reports); the 1e-3 factor in `sed_equilibrium` performs that
mixed-unit conversion.

Relation to the Fortran (sed/src/radfield.f90, sed/src/p_sub.f90):

  * `planck_lambda` is `radfield :: bbody`, with the same physical constants
    (CODATA 2010) and the same x >= 700 cutoff; exp(x) - 1 is np.expm1
    throughout, where bbody switches to a sinh form only below x = 1e-4.
  * `absorbed_power` is `p_sub :: radiation_absorption`: the logarithmic
    trapezoid on the field's significant support
    (`radfield :: radiation_field_support`).
  * `calc_P` builds the same transition matrix as `p_sub :: calc_P` and solves
    the same recursion; given the same kappB the two agree to rounding. The
    kappB used here is `emitted_power` (Simpson's rule on the input grid),
    whereas sed_astrodust_mod :: build_kappB integrates C_abs B_lambda on a
    separate logarithmic Planck grid, so a P(T) from the full Fortran solver
    also carries that quadrature difference.
  * `equilibrium_T` root-finds the continuous `emitted_power`; the Fortran
    calc_Teq interpolates log kappB in log T on its temperature grid.

Why a Python copy: this is the "downstream" form of the dust SED
solver -- the function signature takes only what a 3D-RT cell needs to
provide (`lam`, `J_lam`, `C_abs`) and returns what the cell needs to
consume (`lam * I_lam`). Writing it in Python first makes it easy to
prototype radiation field shapes, plot intermediate quantities, and
validate the Fortran library on representative grains without the
full T-matrix and size-distribution calculation.
"""

from __future__ import annotations
import math
from typing import NamedTuple, Tuple

import numpy as np
from scipy.optimize import brentq
from scipy.integrate import simpson


# ---------------------------------------------------------------------------
# Physical constants: the literals of sed/src/radfield.f90 and p_sub.f90
# ---------------------------------------------------------------------------

_H_SI = 6.62606957e-34     # J s
_C_SI = 2.99792458e8       # m/s
_KB_SI = 1.3806488e-23     # J/K
_HC2_SI = 2.0 * _H_SI * _C_SI**2
_HC_KB_SI = _H_SI * _C_SI / _KB_SI
_HC_ERG_UM = 6.62606957e-27 * 2.99792458e10 * 1.0e4   # h*c in erg um (p_sub :: hc)
_FOURPI = 4.0 * math.pi

# Relative intensity threshold for the outer field-support bounds
# (radfield :: J_REL_FLOOR); see _radiation_support.
J_REL_FLOOR = 1.0e-12

# Absorption checks rejected because the field has a single significant
# sample (zero bandwidth); radfield :: isolated_radiation_support_count.
# Read it before/after a solve to detect an unresolved single-sample field.
isolated_radiation_support_count = 0


# ---------------------------------------------------------------------------
# Planck function
# ---------------------------------------------------------------------------

def planck_lambda(lam_um: np.ndarray, T: float) -> np.ndarray:
    """B_lambda(T) at lam_um, in W m^-3 sr^-1 (SI, per unit wavelength: W m^-2 m^-1 sr^-1).

    Mirrors the Fortran `radfield :: bbody` (SI): x = hc/(lambda k T) >= 700
    gives 0 (the tail below the double range), and exp(x) - 1 is taken
    without cancellation as x -> 0. Multiply by 10 to get CGS
    (erg s^-1 cm^-3 sr^-1 per cm wavelength).
    """
    lam_m = np.asarray(lam_um, dtype=np.float64) * 1.0e-6
    if T <= 0.0:
        return np.zeros(np.shape(lam_m))[()]
    with np.errstate(divide='ignore', over='ignore', invalid='ignore'):
        x = _HC_KB_SI / (T * lam_m)
        ok = (lam_m > 0.0) & (x < 700.0)
        B = np.where(ok, _HC2_SI / lam_m**5 / np.expm1(np.where(ok, x, 1.0)), 0.0)
    return B[()]


# ---------------------------------------------------------------------------
# Field support and absorbed power
# ---------------------------------------------------------------------------

def _radiation_support(lam_um: np.ndarray, J_lam: np.ndarray) -> Tuple[int, int]:
    """Significant continuum endpoints of the field, as 1-based indices.

    radfield :: radiation_field_support with the grid checks of
    p_sub :: radiation_support. The wavelength grid must be positive and
    strictly increasing. Outer bounds are the shortest and longest samples
    with J_lam > J_REL_FLOOR * max(J_lam); interior zero samples remain
    part of the represented continuum. Returns (0, 0) for a field with no
    significant sample, and also for one with a single significant sample,
    which has zero bandwidth (counted in isolated_radiation_support_count).
    """
    global isolated_radiation_support_count
    if lam_um.shape != J_lam.shape:
        raise ValueError('radiation_support: size mismatch')
    if np.any(lam_um <= 0.0):
        raise ValueError('radiation_support: nonpositive wavelength')
    if lam_um.size > 1 and np.any(lam_um[1:] <= lam_um[:-1]):
        raise ValueError('radiation_support: unordered wavelengths')
    if lam_um.size == 0:
        return 0, 0
    peak = J_lam.max()
    if not peak > 0.0:
        return 0, 0
    significant = np.flatnonzero(J_lam > J_REL_FLOOR * peak)
    if significant.size == 0:
        return 0, 0
    first, last = int(significant[0]) + 1, int(significant[-1]) + 1
    if first == last:
        isolated_radiation_support_count += 1
        return 0, 0
    return first, last


def absorbed_power(lam_um: np.ndarray, Cabs: np.ndarray, J_lam: np.ndarray) -> float:
    """Total absorbed power per grain (mixed-unit: cm^2 * SI Blambda * dlambda-um).

    `int Cabs(lam) J_lam(lam) dlam` by the logarithmic trapezoid rule
    (Cabs J_lam lam linear in ln lam between samples) on the field's
    significant support, as `p_sub :: radiation_absorption`; no half
    interval toward an exterior zero sample contributes. Multiply by 4*pi
    to get the total over solid angle.
    """
    lam_um = np.asarray(lam_um, dtype=np.float64)
    Cabs = np.asarray(Cabs, dtype=np.float64)
    J_lam = np.asarray(J_lam, dtype=np.float64)
    if Cabs.shape != lam_um.shape:
        raise ValueError('radiation_absorption: size mismatch')
    first, last = _radiation_support(lam_um, J_lam)
    if first == 0:
        return 0.0
    s = slice(first - 1, last)
    y = Cabs[s] * J_lam[s] * lam_um[s]
    width = np.diff(np.log(lam_um[s]))
    return float(np.sum(0.5 * (y[:-1] + y[1:]) * width))


def emitted_power(lam_um: np.ndarray, Cabs: np.ndarray, T: float) -> float:
    """`int Cabs(lam) B_lam(T) dlam` (mixed-unit: cm^2 * SI Blambda * dlambda-um).

    Simpson's rule on the input grid. Multiply by 4*pi for the total
    radiated luminosity of the grain.
    """
    return simpson(Cabs * planck_lambda(lam_um, T), x=lam_um)


# ---------------------------------------------------------------------------
# Equilibrium temperature
# ---------------------------------------------------------------------------

def equilibrium_T(lam_um: np.ndarray,
                  Cabs: np.ndarray,
                  J_lam: np.ndarray,
                  T_lo: float = 2.7,
                  T_hi: float = 3.0e3) -> float:
    """Solve `int Cabs J_lam dlam = int Cabs B_lam(Teq) dlam` for Teq.

    The absorbed side is `absorbed_power` (the support and interpolant of
    p_sub :: radiation_absorption). The 4*pi cancels on both sides (it's
    the same on absorption and emission for an isotropic grain in LTE).
    brentq bracketing.
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
# Stochastic heating (Guhathakurta & Draine 1989 matrix solver with the
# Draine & Li 2001 upward rates), p_sub.f90 routine by routine
# ---------------------------------------------------------------------------
#
# The routines below keep the 1-based indices of p_sub.f90 (element 0 of each
# array is unused) so that every index expression can be checked against the
# Fortran line by line.

# 1/k, k = 1 .. 15, for the series of _exponential_remainders
_INV_INTEGER = (None,) + tuple(1.0 / k for k in range(1, 16))


def _log_one_plus(x: float) -> float:
    """ln(1 + x) for x >= 0, as ln(u) x / (u - 1) with u = 1 + x, which
    cancels the rounding of u (accurate to a few ulps at any x)."""
    u = 1.0 + x
    if u == 1.0:
        return x
    return math.log(u) * x / (u - 1.0)


def _exponential_remainders(t: float) -> Tuple[float, float]:
    """phi = e^t - 1 - t = sum_(k>=2) t^k/k! and
    psi = t e^t - (e^t - 1) - t^2/2 = sum_(k>=3) (k-1) t^k/k!  for t >= 0.

    Below t = 0.3 both are summed as series, truncated where the next term
    is under 1e-17 of the sum; above it the closed forms lose to
    cancellation at most a factor 40 in relative rounding (psi at 0.3).
    """
    if t < 0.3:
        # highest power kept; psi, whose sum starts at t^3, sets it
        if t < 1.0e-4:
            n = 6
        elif t < 1.0e-3:
            n = 7
        elif t < 1.0e-2:
            n = 9
        elif t < 0.1:
            n = 12
        else:
            n = 15
        term = 0.5 * t * t
        phi = term
        psi = 0.0
        for k in range(3, n + 1):
            term = term * t * _INV_INTEGER[k]
            phi = phi + term
            psi = psi + float(k - 1) * term
        return phi, psi
    e = math.exp(t)
    return e - 1.0 - t, t * e - (e - 1.0) - 0.5 * t * t


class _AbsorbedPowerSpectrum:
    """The absorbed power per unit photon energy on a field's support
    (p_sub :: absorbed_power_t), arrays 1-based of length m."""
    __slots__ = ('m', 'E', 'y', 'slope', 'F0', 'F1')


class _PhotonEnergy(NamedTuple):
    """A photon energy E [erg] on that spectrum: its segment i, y there, and
    the absorbed power F0 = int_0^E g dE below it (p_sub :: photon_energy_t)."""
    E: float
    y: float
    F0: float
    i: int


def _absorbed_power_spectrum(lam_um, Cabs, J_lam, first, last) -> _AbsorbedPowerSpectrum:
    """The absorbed power per unit photon energy on the field support, in
    ascending photon energy E_n = hc/lambda: g(E) dE = 4 pi y d(ln E),
    y = Cabs J_lambda lambda, linear in ln(E) between samples -- the
    interpolant of absorbed_power's logarithmic trapezoid. The segment
    widths tau_n in ln(E) are taken from the sample energies themselves,
    ln(1 + relative step), the scale on which every photon energy inside a
    segment is placed (_photon_energy); the total power then agrees with
    absorbed_power to the rounding of its widths. Stores the slopes
    dy/d(ln E) and the cumulants F0(E_n) = int_0^E_n g dE and
    F1(E_n) = int_0^E_n E g dE.
    """
    g = _AbsorbedPowerSpectrum()
    m = last - first + 1
    g.m = m
    E = [0.0] * (m + 1)
    y = [0.0] * (m + 1)
    slope = [0.0] * (m + 1)
    F0 = [0.0] * (m + 1)
    F1 = [0.0] * (m + 1)
    for n in range(1, m + 1):
        k = last - n           # 0-based index of Fortran's k = last - n + 1
        lam_k = float(lam_um[k])
        E[n] = _HC_ERG_UM / lam_k
        y[n] = float(Cabs[k]) * float(J_lam[k]) * lam_k
    for n in range(1, m):
        tau = _log_one_plus((E[n+1] - E[n]) / E[n])
        slope[n] = (y[n+1] - y[n]) / tau
        phi, psi = _exponential_remainders(tau)
        F0[n+1] = F0[n] + _FOURPI * (y[n] + 0.5 * slope[n] * tau) * tau
        F1[n+1] = F1[n] + _FOURPI * E[n] * (y[n] * (phi + tau) + slope[n] * (psi + 0.5 * tau * tau))
    g.E, g.y, g.slope, g.F0, g.F1 = E, y, slope, F0, F1
    return g


def _photon_energy(g: _AbsorbedPowerSpectrum, x: float, hint: int) -> _PhotonEnergy:
    """Photon energy x [erg] placed on the absorbed-power spectrum: its
    segment i (E(i) <= x < E(i+1); 0 below the softest sample, m at or
    above the hardest), y there, and F0(x). The search steps from hint,
    since successive energies lie close together."""
    E = g.E
    if x < E[1]:
        return _PhotonEnergy(x, 0.0, 0.0, 0)
    if x >= E[g.m]:
        return _PhotonEnergy(x, g.y[g.m], g.F0[g.m], g.m)
    i = min(max(hint, 1), g.m - 1)
    while x < E[i]:
        i -= 1
    while x >= E[i+1]:
        i += 1
    t = 0.0
    if x > E[i]:
        t = _log_one_plus((x - E[i]) / E[i])
    y = g.y[i] + g.slope[i] * t
    F0 = g.F0[i] + _FOURPI * (g.y[i] + 0.5 * g.slope[i] * t) * t
    return _PhotonEnergy(x, y, F0, i)


def _segment_moments(g, i, start, y_start, width) -> Tuple[float, float]:
    """int g dE and int (E - start) g dE over (start, start + width) inside
    segment i, in closed form for y linear in ln(E): with
    t = ln(1 + width/start), y_s = y at start and k the segment slope,
       d0 = 4 pi (y_s + k t/2) t,
       r1 = 4 pi start [y_s (e^t - 1 - t) + k (t e^t - (e^t - 1) - t^2/2)].
    """
    t = _log_one_plus(width / start)
    phi, psi = _exponential_remainders(t)
    d0 = _FOURPI * (y_start + 0.5 * g.slope[i] * t) * t
    r1 = _FOURPI * start * (y_start * phi + g.slope[i] * psi)
    return d0, r1


def _absorption_moments(g, p: _PhotonEnergy, q: _PhotonEnergy, w: float) -> Tuple[float, float]:
    """d0 = int g dE and r1 = int (E - p) g dE over the photon energies
    (p, q), q - p = w >= 0, the part outside the field support contributing
    nothing. The interval is split at the samples it crosses: partial
    segments are integrated in closed form from their start and width,
    whole segments from the cumulants. q only locates the end; the width
    is w, so the accuracy of a narrow interval is that of w, not that of
    the difference of two photon energies."""
    if w <= 0.0 or p.i == g.m or q.i == 0:
        return 0.0, 0.0
    E = g.E
    start = p.E
    y_start = p.y
    i = p.i
    offset = 0.0
    if i == 0:
        # below the softest field photon the integrand is zero
        offset = E[1] - p.E
        start = E[1]
        y_start = g.y[1]
        i = 1
    width = w - offset
    j = q.i
    if j == g.m:
        # above the hardest field photon the integrand is zero
        j = g.m - 1
        width = E[g.m] - start

    if j <= i:
        d0, r1 = _segment_moments(g, i, start, y_start, max(width, 0.0))
    else:
        d0, r1 = _segment_moments(g, i, start, y_start, E[i+1] - start)
        if j > i + 1:
            # whole segments i+1 .. j-1
            piece0 = g.F0[j] - g.F0[i+1]
            r1 = r1 + ((g.F1[j] - g.F1[i+1]) - E[i+1] * piece0) + (E[i+1] - start) * piece0
            d0 = d0 + piece0
        piece0, piece1 = _segment_moments(g, j, E[j], g.y[j], max(width - (E[j] - start), 0.0))
        r1 = r1 + piece1 + (E[j] - start) * piece0
        d0 = d0 + piece0
    r1 = r1 + offset * d0
    return d0, r1


def _upward_absorption_rates(lam_um, Cabs, J_lam, first, last, H) -> np.ndarray:
    """Photon-absorption rates A(u, l) [s^-1] from enthalpy bin l to every
    higher bin u, Draine & Li (2001, ApJ 551, 807; DL01) eqs. (15)-(24), (16)
    and (28); p_sub :: upward_absorption_rates. Returns the (nT, nT) matrix
    with A(u, l) in row u, column l (0-based), zero on and above the
    diagonal.

    Bin u is [E_u,min, E_u,max] = [edge(u-1), edge(u)], with each inner edge
    at the geometric mean sqrt(H_u H_(u+1)) of adjacent enthalpies and the
    two outer edges one half step beyond the end points, so the width dH_u
    agrees with GD89's H_u ln(H_(u+1)/H_(u-1))/2 to second order in the
    step.

    Approximation: the representative energy of bin u is the grid enthalpy
    H_u, at which T(u), kappB(u) and the cooling difference quotient are
    defined, and H_u - H_l is the energy of a transition. DL01 Sec. 4.1
    instead takes E_u = (E_u,min + E_u,max)/2. On a grid given in advance
    these cannot both hold with contiguous bins unless it is uniform; the
    geometric edges center H_u in ln H. With r- = H_u/H_(u-1) and
    r+ = H_(u+1)/H_u, the arithmetic center is
       H_u (r-^(-1/2) + r+^(1/2))/2
         = H_u [1 + (ln r+ - ln r-)/4 + (ln r)^2/8 + ...],
    where the first term vanishes on a uniform ln H grid and is of second
    order where the spacing in ln H varies smoothly, as it does on H(T).
    The relative offset (ln r)^2/8 is 1.2e-5 at r = 1.01, 1.1e-3 at
    r = 1.1, 4.2e-3 at r = 1.2 and 2.1e-2 at r = 1.5. It leaves the energy
    budget exact (absorption and cooling use the same H), shifts each rate
    by a relative amount of that order, and falls quadratically when the
    enthalpy grid is refined; that refinement check (see calc_P) is what
    bounds it for a given window.

    With g(E) dE = 4 pi Cabs J_lambda dlambda, the power absorbed per unit
    photon energy, eq. (15) is
       A(u, l) = [int Phi_ul(E) g(E) dE] / (H_u - H_l),
    where Phi_ul(E) = dH_u G_ul(E) is the fraction of a grain spread
    uniformly over bin l that a photon E moves into bin u (eqs. 17-24):
    (E - W1)/dH_l on (W1, W2), min(dH_l, dH_u)/dH_l on (W2, W3) and
    (W4 - E)/dH_l on (W3, W4), with W1 = E_u,min - E_l,max,
    W2 = W1 + min(dH_l, dH_u), W3 = W1 + max(dH_l, dH_u) and
    W4 = W1 + dH_l + dH_u. In terms of A_k = edge(k) - edge(l) and
    B_k = edge(k) - edge(l-1), W1 = A_(u-1), W4 = B_u and {W2, W3} =
    {B_(u-1), A_u}, so stepping u adds the two photon energies A_u, B_u.

    Eq. (16): the top bin nT takes every photon that carries the grain past
    its lower edge, with weight (E - W1)/dH_l up to Wc = W1 + dH_l and 1
    beyond. Eq. (28): photons E < dH_l leave a fraction 1 - E/dH_l of bin l
    inside bin l; that power is added to l -> l+1. Summed over u these
    fractions are 1 at every E > 0, so each bin l < nT absorbs the whole
    power, 4 pi absorbed_power(lam_um, Cabs, J_lam).

    The rising and falling edges are integrated in closed form from their
    start over their width (_absorption_moments), not as differences of the
    cumulants F0 and F1 at their ends: a bin can be many orders of
    magnitude narrower than the photon energies, and the moment
    int (E - W) g dE taken as F1 - W F0 at E ~ W would lose the digits of
    E/dH. The plateau and the top-bin tail are differences of F0, whose
    rounding is a fraction of order 1e-16 of the bin's whole absorption.
    """
    nT = len(H)
    Amat = np.zeros((nT + 1, nT + 1))
    if first == 0 or last <= first or nT < 2:
        return Amat[1:, 1:]
    g = _absorbed_power_spectrum(lam_um, Cabs, J_lam, first, last)
    Hp = [0.0] + [float(v) for v in H]
    sqrt = math.sqrt

    edge = [0.0] * (nT + 1)
    edge[0] = Hp[1] * sqrt(Hp[1] / Hp[2])
    for l in range(1, nT):
        edge[l] = sqrt(Hp[l] * Hp[l+1])
    edge[nT] = Hp[nT] * sqrt(Hp[nT] / Hp[nT-1])
    # The widths edge(l) - edge(l-1) from the enthalpy differences, which
    # are exact, rather than from the difference of two edges, which loses
    # the digits of H/dH on a narrow window.
    dH = [0.0] * (nT + 1)
    dH[1] = sqrt(Hp[1] / Hp[2]) * (Hp[2] - Hp[1])
    for l in range(2, nT):
        dH[l] = sqrt(Hp[l]) * (Hp[l+1] - Hp[l-1]) / (sqrt(Hp[l+1]) + sqrt(Hp[l-1]))
    dH[nT] = sqrt(Hp[nT] / Hp[nT-1]) * (Hp[nT] - Hp[nT-1])

    zero = _photon_energy(g, 0.0, 1)
    hardest = _photon_energy(g, g.E[g.m], g.m)
    B_lo_i = 1
    for l in range(1, nT):
        dH_l = dH[l]
        A_lo = zero
        B_lo = _photon_energy(g, dH_l, B_lo_i)
        # eq. (28): weight dH_l - E on (0, dH_l)
        d0, r1 = _absorption_moments(g, zero, B_lo, dH_l)
        Amat[l+1, l] = (dH_l * d0 - r1) / (dH_l * (Hp[l+1] - Hp[l]))
        for u in range(l + 1, nT + 1):
            if A_lo.i == g.m:
                break
            if u == nT:
                # eq. (16): weight E - W1 on (W1, Wc), dH_l on (Wc, E_hard),
                # with Wc = W1 + dH_l = B_(nT-1)
                d0, r1 = _absorption_moments(g, A_lo, B_lo, dH_l)
                power = r1 + dH_l * (hardest.F0 - B_lo.F0)
            else:
                dH_u = dH[u]
                rise = min(dH_l, dH_u)
                A_hi = _photon_energy(g, edge[u] - edge[l], A_lo.i)
                B_hi = _photon_energy(g, edge[u] - edge[l-1], B_lo.i)
                power = 0.0
                if B_hi.i > 0:
                    # rising edge (A_(u-1), W2), plateau (W2, W3), falling
                    # edge (W3, B_u)
                    if B_lo.E <= A_hi.E:
                        d0, r1 = _absorption_moments(g, A_lo, B_lo, rise)
                        power = r1 + rise * max(A_hi.F0 - B_lo.F0, 0.0)
                        d0, r1 = _absorption_moments(g, A_hi, B_hi, rise)
                    else:
                        d0, r1 = _absorption_moments(g, A_lo, A_hi, rise)
                        power = r1 + rise * max(B_lo.F0 - A_hi.F0, 0.0)
                        d0, r1 = _absorption_moments(g, B_lo, B_hi, rise)
                    power = power + (rise * d0 - r1)
                A_lo = A_hi
                B_lo = B_hi
            Amat[u, l] = Amat[u, l] + max(power, 0.0) / (dH_l * (Hp[u] - Hp[l]))
        B_lo_i = B_lo.i
    return Amat[1:, 1:]


def _build_kappB(lam_um: np.ndarray,
                 Cabs: np.ndarray,
                 T_grid: np.ndarray) -> np.ndarray:
    """`int Cabs(lam) B_lam(T) dlam` for each T in T_grid.

    Used by calc_P. Returns array of length len(T_grid) in mixed units
    (cm^2 * SI Blambda * dlambda-um).
    """
    return np.array([emitted_power(lam_um, Cabs, T) for T in T_grid])


def _transition_matrix(lam_um, Cabs, J_lam, H, kappB, kappCMB=0.0) -> np.ndarray:
    """GD89 transition matrix A(f, i) [s^-1] (row f = final bin, column
    i = initial bin, 0-based): the continuous cooling
    A(u-1, u) = 4 pi [kappB(u) - kappCMB] / (H(u) - H(u-1)) above the
    diagonal and the DL01 absorption rates (_upward_absorption_rates)
    below it."""
    first, last = _radiation_support(lam_um, J_lam)
    Amat = _upward_absorption_rates(lam_um, Cabs, J_lam, first, last, H)
    nT = len(H)
    for i1 in range(1, nT):
        kappB_tot = kappB[i1] - kappCMB
        Amat[i1-1, i1] = _FOURPI / (H[i1] - H[i1-1]) * kappB_tot
    return Amat


def _steady_state_distribution(Amat: np.ndarray) -> np.ndarray:
    """Steady state P of the GD89 linear recursion for the transition
    matrix Amat (_transition_matrix), normalized to sum 1. With
    B(j, i) = sum_(k>=j) A(k, i), k > i, the balance across the boundary
    below bin j gives P(j) = sum_(i<j) B(j, i) P(i) / A(j-1, j). Without
    absorbed photons the stationary state is the lowest bin, P(1) = 1."""
    nT = Amat.shape[0]
    lower = np.tril(Amat, -1)
    Bmat = np.tril(np.cumsum(lower[::-1, :], axis=0)[::-1, :], -1)
    P = np.zeros(nT)
    if np.sum(Bmat) > 0.0:
        P[0] = 1.0
        for j in range(1, nT):
            if Amat[j-1, j] > 0.0:
                P[j] = np.dot(Bmat[j, :j], P[:j]) / Amat[j-1, j]
            if P[j] > 1e50:
                P[:j+1] /= P[j]      # rescale to avoid overflow
        sP = P.sum()
        if sP > 0.0:
            P /= sP
    else:
        P[0] = 1.0
    return P


def calc_P(lam_um: np.ndarray,
           Cabs: np.ndarray,
           J_lam: np.ndarray,
           T_grid: np.ndarray,
           U_grain: np.ndarray,
           kappCMB: float = 0.0) -> np.ndarray:
    """Temperature distribution P(T) of a stochastically heated grain:
    the Guhathakurta & Draine (1989, GD89) recursion with the finite-bin
    upward rates of Draine & Li (2001, DL01).

    Follows `p_sub :: calc_P` step by step: downward rates are GD89's
    continuous cooling, A(u-1, u) = 4 pi [kappB(u) - kappCMB] /
    (H(u) - H(u-1)); upward rates are DL01 eqs. (15)-(24) for l -> u < nT,
    eq. (16) for the top bin nT, and eq. (28) for the absorptions that stay
    inside bin l, added to l -> l+1 (see _upward_absorption_rates), on the
    field's significant support with C_abs J_lambda lambda linear in
    ln lambda between samples, integrated in closed form. For bins narrow
    compared with the photon energies and with the scale on which C_abs J
    varies, these reduce to DL01 eq. (26), i.e. GD89 eq. (2.5),
    A(u, l) = 4 pi C_abs J_lambda hc Delta H_u / E^3, E = H_u - H_l; that
    point-sampled limit fails wherever a bin is broader than the band the
    field occupies.

    Validity: summed over u, every bin l < nT absorbs the power
    4 pi absorbed_power(lam_um, Cabs, J_lam) to rounding, so the steady
    state emits 4 pi sum_(u>=2) P(u) [kappB(u) - kappCMB] = (1 - P(nT))
    times that power. The top bin has no upward transition, so the
    enthalpy window must leave P(nT) negligible. What finite bins still
    approximate is the uniform spread within a bin, the cooling difference
    quotient and the representative energy H_u of a bin (not DL01's
    arithmetic center; see _upward_absorption_rates), so the emitted
    spectrum's shape must be checked by refining the enthalpy grid. A field
    that absorbs nothing leaves P(1) = 1.

    Given the same kappB, the result equals the Fortran calc_P to rounding;
    kappB here is `emitted_power` on the input grid (see the module
    docstring).

    Parameters
    ----------
    lam_um   : (NW,) wavelength grid [um], positive and strictly increasing
    Cabs     : (NW,) absorption cross section [cm^2] for the grain
    J_lam    : (NW,) mean intensity on the same lam_um grid [SI W/m^3/sr]
    T_grid   : (NT,) temperature grid [K], log-spaced and ideally already
               narrowed around the expected Teq for this grain
    U_grain  : (NT,) enthalpy U(T_grid) [erg], strictly increasing
    kappCMB  : scalar, CMB Planck integral at `radfield :: cmb_temperature()`
               -- 2.725 K unless `use_mathis_corrected` is off, which
               restores Mathis (1983)'s 2.9 K. Subtracted from kappB, so
               kappB - kappCMB is the net cooling rate and a grain cannot
               cool below the CMB the field carries. Pass 0 to skip.

    Returns
    -------
    P : (NT,) probability on T_grid (sums to 1).
    """
    lam_um = np.asarray(lam_um, dtype=np.float64)
    Cabs = np.asarray(Cabs, dtype=np.float64)
    J_lam = np.asarray(J_lam, dtype=np.float64)
    T_grid = np.asarray(T_grid, dtype=np.float64)
    U_grain = np.asarray(U_grain, dtype=np.float64)
    NW = len(lam_um)
    NT = len(T_grid)
    if Cabs.shape != (NW,) or J_lam.shape != (NW,):
        raise ValueError('calc_P: Cabs and J_lam must match lam_um')
    if U_grain.shape != (NT,):
        raise ValueError('calc_P: U_grain must match T_grid')

    kappB = _build_kappB(lam_um, Cabs, T_grid)
    Amat = _transition_matrix(lam_um, Cabs, J_lam, U_grain, kappB, kappCMB)
    return _steady_state_distribution(Amat)


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
        Temperature probability on T_grid (sums to 1).
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

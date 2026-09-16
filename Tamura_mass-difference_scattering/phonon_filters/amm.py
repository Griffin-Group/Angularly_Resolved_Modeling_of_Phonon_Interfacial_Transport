"""
Acoustic Mismatch Model (AMM).

Phonons are treated as plane waves undergoing specular refraction at planar
interfaces (Snell's law) with elastic, energy-conserving transmission set by
the acoustic-impedance contrast (paper Eq. 2). Each interface couples one
polarization (longitudinal 'L' or transverse 'T') characterised by the layer
density rho and the corresponding sound speed c (Z = rho c).

Provided here:

* ``interface_power_transmission`` -- the angle-resolved single-interface
  energy transmission alpha(mu) = 1 - R (drives Figs. 5 and 6).
* ``effective_transmission`` -- the angle-averaged single number
  T_eff = 2 * integral_0^1 alpha(mu) mu dmu (drives Fig. S1).
* ``stack_transmission_isotropic`` -- the multilayer angular transmission used
  to combine acoustic mismatch across the Si/SiGe/Al stack with the mass-
  difference scattering (drives the AMM part of Figs. 2 and 3).

Snell refraction: sin(theta_t) = (c2/c1) sin(theta_i); when this exceeds 1 the
phonon is totally internally reflected (alpha = 0) -- the critical-angle cutoff
visible in Figs. 5b, 5d, and 6.

References
----------
Acoustic-mismatch interface transmission alpha(mu) = 1 - R (manuscript Eq. 2)
and the angle-averaged effective transmission T_eff:

* W. A. Little, "The transport of heat between dissimilar solids at low
  temperatures," Can. J. Phys. 37, 334-349 (1959).
  https://doi.org/10.1139/p59-037  [original acoustic-mismatch model]
* E. T. Swartz and R. O. Pohl, "Thermal boundary resistance," Rev. Mod. Phys.
  61, 605-668 (1989).  https://doi.org/10.1103/RevModPhys.61.605  [AMM review]

General nanoscale-transport treatment: G. Chen, "Nanoscale Energy Transport and
Conversion" (Oxford Univ. Press, 2005), Ch. 5.
"""

from __future__ import annotations

import numpy as np


def _impedance(layer: dict, mode: str) -> float:
    """Acoustic impedance Z = rho * c for the requested polarization."""
    c = layer["cL"] if mode == "L" else layer["cT"]
    return layer["rho"] * c


def _speed(layer: dict, mode: str) -> float:
    return layer["cL"] if mode == "L" else layer["cT"]


def interface_power_transmission(mat1: dict, mat2: dict, mode: str, mu):
    """Single-interface energy transmission alpha(mu) from mat1 into mat2.

    With mu = cos(theta_i) and the Snell-refracted mu_t = cos(theta_t),

        t     = 2 Z1 mu / (Z2 mu_t + Z1 mu)        (pressure transmission)
        alpha = t^2 (Z2 mu_t) / (Z1 mu)            (= 1 - R, Eq. 2)

    Totally internally reflected angles (sin theta_t >= 1) return alpha = 0.

    Parameters
    ----------
    mat1, mat2 : dict
        Layer property dicts with keys 'rho', 'cL', 'cT'.
    mode : {'L', 'T'}
        Phonon polarization.
    mu : float or ndarray
        Cosine of the incidence angle in ``mat1``.

    Returns
    -------
    ndarray
        Energy transmission probability alpha in [0, 1].
    """
    mu = np.atleast_1d(np.asarray(mu, dtype=float))
    c1, c2 = _speed(mat1, mode), _speed(mat2, mode)
    Z1, Z2 = _impedance(mat1, mode), _impedance(mat2, mode)

    theta_i = np.arccos(np.clip(mu, -1.0, 1.0))
    sin_t = (c2 / c1) * np.sin(theta_i)

    alpha = np.zeros_like(mu)
    # Transmitted where not totally internally reflected (sin_t < 1) and the
    # incidence is not exactly grazing (mu > 0; at mu = 0 there is no normal
    # energy flux, so alpha = 0 and we avoid a 0/0 in the power ratio).
    transmitted = (sin_t < 1.0) & (mu > 0.0)
    mu_t = np.cos(np.arcsin(sin_t[transmitted]))
    m = mu[transmitted]
    t = (2 * Z1 * m) / (Z2 * mu_t + Z1 * m)
    alpha[transmitted] = t**2 * (Z2 * mu_t) / (Z1 * m)
    return alpha


def effective_transmission(mat1: dict, mat2: dict, mode: str, n_mu: int = 4000) -> float:
    """Angle-averaged effective transmission across a single interface (Fig. S1).

    T_eff = 2 * integral_0^1 alpha(mu) mu dmu

    weights each incidence angle by its contribution (mu) to energy transfer,
    yielding one frequency-independent number per interface and polarization.

    Parameters
    ----------
    mat1, mat2 : dict
        Layer property dicts.
    mode : {'L', 'T'}
        Polarization.
    n_mu : int
        Number of mu quadrature points on [0, 1].

    Returns
    -------
    float
        Effective transmission T_eff.
    """
    mu = np.linspace(0.0, 1.0, n_mu)
    alpha = interface_power_transmission(mat1, mat2, mode, mu)
    return float(2.0 * np.trapezoid(alpha * mu, mu))


def stack_transmission_isotropic(angles_rad, stack, mode="L"):
    """Angle-resolved transmission through a multilayer stack (isotropic model).

    Intermediate interfaces are angle-integrated and folded into a scalar
    (capturing the SI's uniform-redistribution boundary model, S5), while the
    final interface retains its full mu dependence. The leading factor 1/2 is
    the incident-side normalisation. This produces the angle-dependent acoustic
    transmission that is later combined with the mass-difference scattering
    spectrum to form Figs. 2 and 3.

    Parameters
    ----------
    angles_rad : ndarray
        Incident angles (radians).
    stack : list of dict
        Ordered layer property dicts (propagation direction).
    mode : {'L', 'T'}
        Polarization.

    Returns
    -------
    T_final_mu : ndarray
        Transmission versus mu (ascending mu), for the full stack.
    mus : ndarray
        Corresponding mu = cos(theta) grid (near-grazing mu <= 1e-5 dropped).
    """
    mus = np.cos(np.asarray(angles_rad))
    mus = np.sort(mus)
    mus = mus[mus > 1e-5]  # drop near-grazing incidence

    speeds = [_speed(layer, mode) for layer in stack]
    impedances = [_impedance(layer, mode) for layer in stack]

    cumulative_scalar = 1.0 / 2.0
    T_final_mu = None

    for i in range(len(stack) - 1):
        Z1, Z2 = impedances[i], impedances[i + 1]
        c1, c2 = speeds[i], speeds[i + 1]

        T_mu = np.empty_like(mus)
        for j, mu in enumerate(mus):
            theta_i = np.arccos(mu)
            sin_t = (c2 / c1) * np.sin(theta_i)
            if abs(sin_t) >= 1.0:
                T_mu[j] = 0.0
                continue
            mu_t = np.cos(np.arcsin(sin_t))
            t = (2 * Z1 * mu) / (Z2 * mu_t + Z1 * mu)
            T_mu[j] = t**2 * (Z2 * mu_t) / (Z1 * mu)

        if len(stack) == 2:
            # Single interface: return the bare angle-resolved transmission.
            return T_mu[::-1], mus

        if i < len(stack) - 2:
            # Intermediate interface: angle-integrate into the running scalar.
            cumulative_scalar *= 2 * np.trapezoid(T_mu, mus)
        else:
            # Final interface: keep full mu dependence.
            T_final_mu = cumulative_scalar * T_mu

    return T_final_mu[::-1], mus

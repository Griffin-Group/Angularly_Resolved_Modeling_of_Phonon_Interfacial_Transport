"""
Thermal-occupation quantities.

  * Bose-Einstein occupation n(omega, T) and its temperature derivative,
  * the density-of-states weighted occupation DOS * n (the grey overlay in
    Fig. 2b, indicating which phonons are thermally populated at a given T),
  * the dimensionless Debye specific-heat integrand and its cumulative
    integral (Fig. S5).

References
----------
The Bose-Einstein occupation and the Debye specific-heat integrand are standard
statistical-mechanics / lattice-dynamics results; see e.g. N. W. Ashcroft and
N. D. Mermin, "Solid State Physics" (1976), Ch. 23, or G. Chen, "Nanoscale
Energy Transport and Conversion" (Oxford Univ. Press, 2005).
"""

from __future__ import annotations

import numpy as np
from scipy import constants as _const

from . import dispersion

HBAR = _const.hbar
K_B = _const.k


def bose_einstein(omega, T: float):
    """Bose-Einstein occupation and its temperature derivative.

    Parameters
    ----------
    omega : float or ndarray
        Angular frequency (rad/s).
    T : float
        Temperature (K), must be > 0.

    Returns
    -------
    n : ndarray
        Occupation number  n = 1 / (exp(x) - 1),  x = hbar*omega / (k_B T).
    dn_dT : ndarray
        Temperature derivative  dn/dT = (x/T) exp(x) / (exp(x) - 1)^2.
    """
    if T <= 0:
        raise ValueError("Temperature must be greater than zero.")
    omega = np.asarray(omega, dtype=float)
    x = np.clip(HBAR * omega / (K_B * T), -599, 599)  # clip guards against overflow
    exp_x = np.exp(x)
    expm1_x = np.clip(np.expm1(x), -1e150, 1e150)
    n = 1.0 / expm1_x
    dn_dT = (x / T) * exp_x / (expm1_x * expm1_x)
    return n, dn_dT


def dos_times_occupation(element: str, omega, T: float, elemental_data: dict):
    """Density of states and its Bose-Einstein weighting for one material.

    Parameters
    ----------
    element : str
        Element name (key into ``elemental_data``).
    omega : ndarray
        Angular frequency grid (rad/s).
    T : float
        Temperature (K).
    elemental_data : dict
        Per-element property table (see ``materials.load_elemental_properties``).

    Returns
    -------
    dos : ndarray
        BvK density of states.
    n : ndarray
        Bose-Einstein occupation.
    dos_n : ndarray
        ``dos * n`` (thermally weighted DOS; the Fig. 2b overlay).
    dos_dn_dT : ndarray
        ``dos * dn/dT`` (spectral heat-capacity weighting).
    """
    atomic_mass, density, v_s, atoms_per_cell, *_ = elemental_data[element]
    N = dispersion.number_density(density, atomic_mass, atoms_per_cell)
    k_0 = dispersion.zone_boundary_wavevector(N)
    omega_0 = dispersion.bvk_cutoff_frequency(k_0, v_s)

    dos = dispersion.density_of_states(omega, omega_0, N)
    n, dn_dT = bose_einstein(omega, T)
    return dos, n, dos * n, dos * dn_dT


def debye_integrand_cumulative(x_max: float = 15.0, n_points: int = 5000):
    """Debye specific-heat integrand and its cumulative integral (Fig. S5).

    The dimensionless integrand is

        f(x) = x^4 e^x / (e^x - 1)^2,     x = hbar*omega / (k_B T),

    whose running integral F(x) = \\int_0^x f(t) dt approaches the Debye
    limit 4 pi^4 / 15 ~= 25.976 as x -> infinity. The cumulative curve shows
    which (dimensionless) phonon frequencies carry the heat capacity.

    Parameters
    ----------
    x_max : float
        Upper limit of the dimensionless frequency grid.
    n_points : int
        Number of grid points.

    Returns
    -------
    x : ndarray
        Dimensionless frequency grid.
    integrand : ndarray
        f(x).
    cumulative : ndarray
        Running integral F(x) (trapezoidal, F(0)=0).
    """
    from scipy.integrate import cumulative_trapezoid

    x = np.linspace(1e-8, x_max, n_points)
    integrand = x**4 * np.exp(x) / (np.exp(x) - 1) ** 2
    cumulative = cumulative_trapezoid(integrand, x, initial=0)
    return x, integrand, cumulative

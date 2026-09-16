"""
Born--von Karman (BvK) phonon dispersion on an isotropic Debye sphere.

The BvK branch is a sine dispersion (paper Eqs. 10-14):

    omega(k) = omega_0 * sin( (pi/2) * k / k_0 )

from which the group velocity and density of states follow analytically. These
quantities feed the Tamura mass-difference scattering model (``scattering.py``)
for the composition-dependent alloy mean free paths (Figs. S4, S6) and the
density-of-states weighting in the thermal-occupation overlays (Fig. 2b).

Symbols
-------
omega   angular frequency           [rad/s]
k       phonon wavevector magnitude [1/m]
k_0     zone-boundary wavevector    [1/m]   k_0 = (6 pi^2 N)^(1/3)
omega_0 BvK cutoff angular freq     [rad/s] omega_0 = (2/pi) v_s k_0
N       number density of primitive cells [1/m^3]
v_s     (averaged) sound velocity   [m/s]

References
----------
The Born-von Karman sine dispersion and the analytic density of states
(manuscript Eqs. 10-14) follow the classical lattice-dynamics model of

* M. Born and Th. von Karman, "Uber Schwingungen in Raumgittern," Physikalische
  Zeitschrift 13, 297-309 (1912).  [no DOI; pre-DOI era]

For a modern textbook treatment see G. Chen, "Nanoscale Energy Transport and
Conversion" (Oxford Univ. Press, 2005), Ch. 3, or N. W. Ashcroft and N. D.
Mermin, "Solid State Physics" (1976), Ch. 22.
"""

from __future__ import annotations

import numpy as np

# atomic mass unit in grams (densities are supplied in g/cm^3)
_AMU_G = 1.6605e-24


def number_density(density_g_cc: float, atomic_mass_amu: float, atoms_per_cell: int) -> float:
    """Number density of primitive unit cells [1/m^3].

    N = rho / (m * amu * n_atoms_per_cell), with the 1e6 factor converting
    cm^-3 to m^-3.
    """
    return density_g_cc / (atomic_mass_amu * _AMU_G * atoms_per_cell) * 1e6


def zone_boundary_wavevector(N: float) -> float:
    """Debye-sphere zone-boundary wavevector k_0 = (6 pi^2 N)^(1/3) [1/m]."""
    return np.cbrt(6 * np.pi**2 * N)


def bvk_cutoff_frequency(k_0: float, v_s: float) -> float:
    """BvK cutoff angular frequency omega_0 = (2/pi) v_s k_0 [rad/s]."""
    return 2.0 * v_s * k_0 / np.pi


def group_velocity(omega, omega_0: float, k_0: float):
    """Group velocity v_g(omega) of the BvK branch [m/s].

    Differentiating omega(k) = omega_0 sin((pi/2) k/k_0):

        v_g = domega/dk = (pi/2)(omega_0/k_0) cos(arcsin(omega/omega_0)).

    Parameters
    ----------
    omega : float or ndarray
        Angular frequency (rad/s); must not exceed ``omega_0``.
    omega_0 : float
        BvK cutoff angular frequency (rad/s).
    k_0 : float
        Zone-boundary wavevector (1/m).

    Raises
    ------
    ValueError
        If any ``omega`` exceeds ``omega_0`` (outside the BvK branch).
    """
    omega = np.asarray(omega, dtype=float)
    if np.any(omega > omega_0):
        # Corrected from the original, which *returned* (rather than raised)
        # the exception object, silently producing a non-array result.
        raise ValueError("All frequencies must be <= the cutoff frequency omega_0.")
    return (np.pi / 2) * (omega_0 / k_0) * np.cos(np.arcsin(omega / omega_0))


def density_of_states(omega, omega_0, N):
    """Phonon density of states D(omega) of the BvK branch (paper Eq. 14).

        D(omega) = (4 k_0^3 / pi^5) * arcsin(omega/omega_0)^2
                   / ( omega_0 * sqrt(1 - (omega/omega_0)^2) )

    Parameters
    ----------
    omega : float or ndarray
        Angular frequency (rad/s).
    omega_0 : float or ndarray
        BvK cutoff angular frequency (rad/s).
    N : float
        Primitive-cell number density (1/m^3), used to set k_0.

    Returns
    -------
    ndarray
        Density of states (states per rad/s per m^3).
    """
    omega = np.asarray(omega, dtype=float)
    omega_0 = np.asarray(omega_0, dtype=float)
    k_0 = np.cbrt(6 * np.pi**2 * N)
    prefactor = 4 * k_0**3 / np.pi**5
    return (
        prefactor
        * np.arcsin(omega / omega_0) ** 2
        / (omega_0 * np.sqrt(1 - (omega / omega_0) ** 2))
    )

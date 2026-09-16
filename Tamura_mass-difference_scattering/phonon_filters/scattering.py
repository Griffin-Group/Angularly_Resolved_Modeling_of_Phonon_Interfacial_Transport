"""
Tamura mass-difference (alloy) phonon scattering.

Random substitution of atoms with different masses on the lattice scatters
phonons at a rate proportional to omega^2 * DOS (Rayleigh-like; paper Eqs. 4-6):

    1/tau(omega) = (pi/6) * V0 * g * omega^p * D(omega)                   (Eq. 4)
    Lambda(omega) = v(omega) / [ (pi/6) V0 g omega^p D(omega) ]          (Eq. 5)
    g = sum_i x_i * ((m_i - m_bar) / m_bar)^2                            (Eq. 6)

Two complementary routes appear in the manuscript and both live here:

* ``mass_difference_mfp_dft`` -- uses the *DFT-computed* density of states
  D(omega) directly, with the physically grounded exponent p = 2 (the Tamura
  rate is intrinsically proportional to omega^2 * D). This drives the SiGe
  results in Figs. 2-4, S7, S8.

* ``alloy_lifetime`` / ``alloy_mfp_min_composition`` / ``mfp_vs_composition``
  -- use the *analytic BvK* density of states, with the empirically fitted
  exponent p = n_aly (2.047 for the SiGe MFP-vs-composition curve in Fig. S4,
  2.00 for the broad alloy survey in Fig. S6). The fitted exponent substitutes
  for omega^2 only on the analytic-DOS branch.

The two exponents are a deliberate modelling choice (real DOS -> p=2; analytic
DOS -> fitted p), not an inconsistency.

References
----------
Mass-difference (alloy/isotope) phonon-scattering rate and mass-fluctuation
factor g (manuscript Eqs. 4-6):

* S. Tamura, "Isotope scattering of dispersive phonons in Ge," Phys. Rev. B 27,
  858-866 (1983).  https://doi.org/10.1103/PhysRevB.27.858

The omega^4 Rayleigh form for point-defect scattering originates with
P. G. Klemens, "The scattering of low-frequency lattice waves by static
imperfections," Proc. Phys. Soc. A 68, 1113-1128 (1955).
https://doi.org/10.1088/0370-1298/68/12/303
"""

from __future__ import annotations

import numpy as np

from . import dispersion
from .dataio import round_sig

# atomic mass unit in kg (per-element masses below are supplied in amu)
_AMU_KG = 1.66e-27


def mass_fluctuation_factor(masses, fractions) -> float:
    """Mass-fluctuation factor g = sum_i x_i ((m_i - m_bar)/m_bar)^2 (Eq. 6).

    Parameters
    ----------
    masses : sequence of float
        Atomic masses of the constituents (any consistent unit).
    fractions : sequence of float
        Atomic fractions x_i (should sum to 1).
    """
    masses = np.asarray(masses, dtype=float)
    fractions = np.asarray(fractions, dtype=float)
    m_bar = np.sum(fractions * masses)
    return float(np.sum(fractions * ((masses - m_bar) / m_bar) ** 2))


def mass_difference_mfp_dft(
    freq_THz,
    dos_per_cell,
    m_Si: float = 28.08 * _AMU_KG,
    m_Ge: float = 72.08 * _AMU_KG,
    V0: float = 48e-30,
    v_La: float = 6281.0,
    v_Ta: float = 2937.0,
):
    """Mass-difference scattering rate and MFP from a DFT density of states.

    For a 50/50 SiGe solid solution, using the DFT DOS directly:

        1/tau = (pi/6) * omega^2 * g * D_cell * (2*pi*1e12)
        Lambda = v_avg / (1/tau)

    where omega is in THz (so ``freq_THz**2`` is the omega^2 factor up to the
    angular conversion folded into the trailing 2*pi*1e12), ``v_avg`` is the
    polarization-averaged sound velocity, and ``g`` follows Eq. 6.

    Parameters
    ----------
    freq_THz : ndarray
        Frequency grid in THz.
    dos_per_cell : ndarray
        Phonon density of states per unit cell on ``freq_THz``.
    m_Si, m_Ge : float
        Atomic masses (kg) of the two species.
    V0 : float
        Volume per unit cell (m^3). Retained for interface compatibility; the
        rate prefactor below follows the original DFT-path normalisation.
    v_La, v_Ta : float
        Longitudinal and transverse sound velocities (m/s).

    Returns
    -------
    inv_tau : ndarray
        Scattering rate (1/s), floored at 1e-20 to avoid division blow-ups.
    mfp : ndarray
        Mean free path (m).
    """
    freq_THz = np.asarray(freq_THz, dtype=float)

    # Eq. 6 mass-fluctuation factor for the 50/50 alloy.
    g = mass_fluctuation_factor([m_Si, m_Ge], [0.5, 0.5])

    # Polarization-averaged sound velocity: 1/v^3 = (1/3)(2/vT^3 + 1/vL^3).
    v_avg = ((1.0 / 3.0) * (2.0 / v_Ta**3 + 1.0 / v_La**3)) ** (-1.0 / 3.0)

    # Tamura rate with the DFT DOS and intrinsic omega^2 dependence.
    inv_tau = (np.pi / 6.0) * (freq_THz**2) * g * dos_per_cell * 2 * np.pi * 1e12
    inv_tau = np.maximum(inv_tau, 1e-20)

    mfp = np.maximum(v_avg / inv_tau, 0.0)
    return inv_tau, mfp


def alloy_lifetime(omega, vs_a, vs_b, omega0_a, omega0_b, N_a, N_b, n_aly, x_a, m_a, x_b, m_b):
    """Tamura alloy-scattering lifetime on the analytic BvK branch (Eq. 4).

    Uses composition-averaged primitive-cell properties (mass, number density,
    sound velocity) to build an averaged BvK density of states D_bar(omega),
    then

        1/tau = (pi/6) * V_bar * g * omega^n_aly * D_bar(omega),
        V_bar = 1 / N_bar.

    Parameters
    ----------
    omega : float or ndarray
        Angular frequency (rad/s).
    vs_a, vs_b : float
        Sound velocities of the two constituents (m/s).
    omega0_a, omega0_b : float
        BvK cutoff angular frequencies (rad/s) (unused directly; the averaged
        cutoff is rebuilt from the averaged number density, matching the model).
    N_a, N_b : float
        Primitive-cell number densities (1/m^3).
    n_aly : float
        Frequency exponent (fitted: 2.047 or 2.00 depending on the figure).
    x_a, m_a, x_b, m_b : float
        Atomic fractions and masses of the two constituents.

    Returns
    -------
    tau : ndarray
        Phonon lifetime (s).
    """
    m_bar = x_a * m_a + x_b * m_b
    N_bar = x_a * N_a + x_b * N_b
    vs_bar = x_a * vs_a + x_b * vs_b

    kc_bar = (6 * np.pi**2 * N_bar) ** (1.0 / 3.0)
    omega0_bar = 2.0 / np.pi * vs_bar * kc_bar

    g = x_a * (1 - m_a / m_bar) ** 2 + x_b * (1 - m_b / m_bar) ** 2

    V_bar = 1.0 / N_bar
    dos_bar = dispersion.density_of_states(omega, omega0_bar, N_bar)

    inv_tau = np.pi / 6.0 * V_bar * g * np.power(omega, n_aly) * dos_bar
    return 1.0 / inv_tau


def _element_bvk_params(element, elemental_data):
    """Return (mass_amu, v_s, N, k_0, omega_0) for an element."""
    atomic_mass, density, v_s, atoms_per_cell, *_ = elemental_data[element]
    N = dispersion.number_density(density, atomic_mass, atoms_per_cell)
    k_0 = dispersion.zone_boundary_wavevector(N)
    omega_0 = dispersion.bvk_cutoff_frequency(k_0, v_s)
    return atomic_mass, v_s, N, k_0, omega_0


def alloy_mfp_min_composition(element_a, element_b, elemental_data, n_aly=2.047):
    """Frequency-resolved alloy MFP at the composition that minimises it.

    The minimising composition is found at 100 GHz (the reference frequency in
    the manuscript): scanning x_a over [0.01, 0.99], the alloy MFP at 100 GHz
    is evaluated and its argmin selected. The full omega-dependent lifetime and
    MFP are then returned at that composition. Drives Fig. S6.

    Parameters
    ----------
    element_a, element_b : str
        Element names (keys into ``elemental_data``).
    elemental_data : dict
        Per-element property table.
    n_aly : float
        Frequency exponent (2.00 for the Fig. S6 survey).

    Returns
    -------
    omega : ndarray
        Angular-frequency grid (rad/s).
    tau : ndarray
        Lifetime at the minimising composition (s).
    mfp : ndarray
        Mean free path at the minimising composition (m).
    min_x_a, min_x_b : float
        Minimising atomic fractions (x_b rounded to 2 sig figs, as in the
        original, with x_a = 1 - x_b).
    """
    m_a, vs_a, N_a, k0_a, w0_a = _element_bvk_params(element_a, elemental_data)
    m_b, vs_b, N_b, k0_b, w0_b = _element_bvk_params(element_b, elemental_data)

    # Frequency grid up to the lower of the two cutoffs (with a small margin so
    # the group velocity / DOS arcsin arguments stay < 1).
    delta = 2 * np.pi  # 1 Hz margin in angular units
    w_top = min(w0_a, w0_b) - delta
    omega = np.logspace(np.log10(100.0), np.log10(w_top), 5000)

    # Composition scan to minimise the 100 GHz MFP.
    x_a_scan = np.linspace(0.01, 0.99, 1000)
    omega_100 = 0.1e12 * 2 * np.pi
    mfp_100 = np.full_like(x_a_scan, np.nan)
    if omega_100 < w0_a and omega_100 < w0_b:
        for j, x_a in enumerate(x_a_scan):
            tau_j = alloy_lifetime(omega_100, vs_a, vs_b, w0_a, w0_b, N_a, N_b,
                                   n_aly, x_a, m_a, 1 - x_a, m_b)
            vg_mix = (x_a * dispersion.group_velocity(omega_100, w0_a, k0_a)
                      + (1 - x_a) * dispersion.group_velocity(omega_100, w0_b, k0_b))
            mfp_100[j] = tau_j * vg_mix

    min_index = int(np.nanargmin(mfp_100))
    min_x_b = round_sig(1 - x_a_scan[min_index], 2)
    min_x_a = 1 - min_x_b

    tau = alloy_lifetime(omega, vs_a, vs_b, w0_a, w0_b, N_a, N_b,
                         n_aly, min_x_a, m_a, min_x_b, m_b)
    v_w = (min_x_a * dispersion.group_velocity(omega, w0_a, k0_a)
           + min_x_b * dispersion.group_velocity(omega, w0_b, k0_b))
    mfp = tau * v_w
    return omega, tau, mfp, min_x_a, min_x_b


def mfp_vs_composition(element_a, element_b, elemental_data,
                       freqs_THz=(0.1, 0.7, 1.1), n_aly=2.047):
    """Alloy MFP versus composition at several fixed frequencies (Fig. S4).

    Parameters
    ----------
    element_a, element_b : str
        Element names. The x-axis is the atomic fraction of ``element_b``.
    elemental_data : dict
        Per-element property table.
    freqs_THz : sequence of float
        Frequencies (THz) at which to plot MFP vs composition.
    n_aly : float
        Frequency exponent (2.047 in the manuscript).

    Returns
    -------
    x_b : ndarray
        Atomic fraction of element_b (the plot abscissa).
    mfp_um : list of ndarray
        MFP (micrometres) versus composition, one array per frequency.
    min_x_a, min_x_b : float
        Composition minimising the MFP at the middle frequency (reference).
    """
    m_a, vs_a, N_a, k0_a, w0_a = _element_bvk_params(element_a, elemental_data)
    m_b, vs_b, N_b, k0_b, w0_b = _element_bvk_params(element_b, elemental_data)

    x_a = np.linspace(0.01, 0.99, 1000)
    x_b = 1 - x_a

    mfp_um = []
    for f_THz in freqs_THz:
        omega = f_THz * 1e12 * 2 * np.pi
        tau = alloy_lifetime(omega, vs_a, vs_b, w0_a, w0_b, N_a, N_b,
                             n_aly, x_a, m_a, 1 - x_a, m_b)
        vg_mix = (x_a * dispersion.group_velocity(omega, w0_a, k0_a)
                  + (1 - x_a) * dispersion.group_velocity(omega, w0_b, k0_b))
        mfp_um.append(tau * vg_mix * 1e6)

    # Reference minimum at the middle frequency (index 1), matching the original.
    ref = 1
    min_index = int(np.argmin(mfp_um[ref]))
    min_x_b = x_b[min_index]
    min_x_a = 1 - min_x_b
    return x_b, mfp_um, min_x_a, min_x_b

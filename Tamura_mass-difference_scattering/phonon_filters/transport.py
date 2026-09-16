"""
Phonon radiative transport through the scattering layer (EPRT via DISORT).

The mass-difference scattering layer is treated with the phonon equivalent of
the radiative transport equation (paper Eqs. 7-9). Because the phonon EPRT and
the optical RTE coincide once thermal emission is dropped, the established
discrete-ordinates solver PythonicDISORT (Ho, ref. [15]) is used to obtain the
diffuse reflectance F_u, diffuse transmittance F_d, and direct transmittance
F_s of a slab of optical depth tau_p, single-scattering albedo sigma, and
Henyey-Greenstein phase function with asymmetry g.

The optical depth is the physical thickness over the frequency-dependent mean
free path, tau_p(nu) = L / Lambda(nu) (paper Eq. 1), so a spectrum of tau maps
to a spectrum of transmittance.

This module wraps PythonicDISORT into the three sweeps the manuscript needs:
  * ``flux_albedo_sweep``   -> transmittance/reflectance vs frequency for a set
    of single-scattering albedos (Figs. 4 and S8),
  * ``flux_vs_thickness``   -> component fluxes vs optical depth (Fig. S7),
  * ``alloy_transmittance_database`` -> transmittance vs frequency for many
    alloys (Fig. S6).

References
----------
Equation of phonon radiative transport (EPRT, manuscript Eqs. 7-9) and the
optical depth tau_p = L/Lambda(nu) (Eq. 1):

* A. Majumdar, "Microscale heat conduction in dielectric thin films," J. Heat
  Transfer 115, 7-16 (1993).  https://doi.org/10.1115/1.2910673
* R. Prasher, "Generalized equation of phonon radiative transport," Appl. Phys.
  Lett. 83, 48-50 (2003).  https://doi.org/10.1063/1.1590421

Discrete-ordinates radiative-transfer solver (manuscript ref. [15]):

* D. J. X. Ho, "PythonicDISORT: A Python reimplementation of the Discrete
  Ordinate Radiative Transfer package DISORT," J. Open Source Softw. 9(103),
  6442 (2024).  https://doi.org/10.21105/joss.06442

Henyey-Greenstein scattering phase function:

* L. G. Henyey and J. L. Greenstein, "Diffuse radiation in the galaxy,"
  Astrophys. J. 93, 70-83 (1941).  https://doi.org/10.1086/144246
"""

from __future__ import annotations

import numpy as np
import PythonicDISORT

from . import dataio, scattering

# ---- DISORT numerical defaults (fixed throughout the manuscript) ----
N_QUAD = 16            # number of discrete-ordinate streams
N_MOMENTS = 32         # Legendre moments supplied for the phase function
ASYMMETRY_G = 0.001    # Henyey-Greenstein asymmetry (near-isotropic)
I0 = np.pi             # incident beam intensity F_i = mu0 * I0
PHI0 = np.pi - 0.01    # incident azimuth (no effect on fluxes)


def henyey_greenstein_legendre(g: float = ASYMMETRY_G,
                               n_moments: int = N_MOMENTS,
                               n_layers: int = 100) -> np.ndarray:
    """Legendre moments g^l of the Henyey-Greenstein phase function.

    Returns an ``(n_layers, n_moments)`` array (every layer identical) as
    PythonicDISORT expects per-layer coefficients.
    """
    return np.tile(g ** np.arange(n_moments), (n_layers, 1))


def _slab_fluxes(tau_i, ssa, mu0, leg_coeffs_row, only_flux=True):
    """Diffuse-up, diffuse-down, direct-down fluxes for a single slab.

    Returns the three fluxes already normalised by the incident flux I0*mu0:
    (reflectance, diffuse transmittance, direct transmittance).
    """
    flux_up, flux_down = PythonicDISORT.pydisort(
        tau_i, ssa, N_QUAD, leg_coeffs_row, mu0, I0, PHI0, only_flux=only_flux
    )[1:3]
    refl = flux_up(0) / (I0 * mu0)
    trans_diff = flux_down(tau_i)[0] / (I0 * mu0)
    trans_direct = flux_down(tau_i)[1] / (I0 * mu0)
    return refl, trans_diff, trans_direct


def flux_albedo_sweep(mfp_trunc, albedos, mu0=1.0, L_um=1.0, leg_coeffs=None):
    """Reflectance and transmittance spectra over a set of albedos (Figs. 4, S8).

    For a fixed scattering-layer thickness ``L_um`` the optical depth spectrum
    is tau = L / Lambda(nu); DISORT is solved at each frequency for every
    single-scattering albedo.

    Parameters
    ----------
    mfp_trunc : ndarray
        Mean free path Lambda(nu) (m) on the (truncated) frequency grid.
    albedos : sequence of float
        Single-scattering albedos sigma to sweep.
    mu0 : float
        Cosine of incidence angle (1.0 = normal). Default 1.0.
    L_um : float
        Scattering-layer thickness in micrometres. Default 1.0.
    leg_coeffs : ndarray, optional
        Phase-function Legendre coefficients; default Henyey-Greenstein.

    Returns
    -------
    dict
        ``{albedo: {"transmittance": ndarray, "reflectance": ndarray}}`` with
        reflectance and (diffuse+direct) transmittance normalised to [0, 1].
    """
    tau = L_um * 1e-6 / mfp_trunc
    if leg_coeffs is None:
        leg_coeffs = henyey_greenstein_legendre(n_layers=max(len(tau), 1))

    out = {}
    for ssa in albedos:
        refl = np.zeros_like(tau)
        tdiff = np.zeros_like(tau)
        tdir = np.zeros_like(tau)
        for i in range(len(tau)):
            refl[i], tdiff[i], tdir[i] = _slab_fluxes(tau[i], ssa, mu0, leg_coeffs[i, :])
        out[ssa] = {
            "transmittance": np.real(tdiff + tdir),
            "reflectance": np.real(refl),
        }
    return out


def transmittance_matrix(mfp_trunc, mu0_array, L_um=10.0, ssa=1 - 1e-6, leg_coeffs=None):
    """Mass-difference transmittance over (angle, frequency) (Fig. 2 input).

    Builds the scattering-layer transmittance spectrum at each incidence angle,
    i.e. the frequency-dependent fraction of phonons surviving the mass-
    difference scattering in the layer, for a near-elastic albedo. This grid is
    later multiplied by the angle-dependent AMM interface transmission to form
    the cumulative transmittance of Fig. 2.

    Parameters
    ----------
    mfp_trunc : ndarray
        Mean free path Lambda(nu) (m) on the truncated frequency grid.
    mu0_array : ndarray
        Incidence cosines (angles) to evaluate.
    L_um : float
        Scattering-layer thickness (micrometres). Default 10.
    ssa : float
        Single-scattering albedo. Default 1 - 1e-6 (near-elastic), matching the
        original Fig. 2 construction.
    leg_coeffs : ndarray, optional
        Phase-function coefficients; default Henyey-Greenstein.

    Returns
    -------
    ndarray, shape (len(mu0_array), len(mfp_trunc))
        Transmittance per (angle, frequency); rows follow ``mu0_array`` order.
    """
    tau = L_um * 1e-6 / mfp_trunc
    if leg_coeffs is None:
        leg_coeffs = henyey_greenstein_legendre(n_layers=max(len(tau), 1))
    matrix = np.zeros((len(mu0_array), len(tau)))
    for k, mu0 in enumerate(mu0_array):
        for i in range(len(tau)):
            _, tdiff, tdir = _slab_fluxes(tau[i], ssa, mu0, leg_coeffs[i, :])
            matrix[k, i] = np.real(tdiff + tdir)
    return matrix


def flux_vs_thickness(tau, mu0_values=(1.0, 0.1), ssa=0.999999, leg_coeffs=None):
    """Component fluxes versus optical depth (Fig. S7).

    Parameters
    ----------
    tau : ndarray
        Optical-depth (acoustic thickness) grid tau_p.
    mu0_values : sequence of float
        Incidence cosines to evaluate (default 1.0 and 0.1).
    ssa : float
        Single-scattering albedo (default 0.999999, near-elastic).
    leg_coeffs : ndarray, optional
        Phase-function coefficients; default Henyey-Greenstein.

    Returns
    -------
    dict
        For each mu0: ``{"F_u","F_d","F_s"}`` (component fluxes, scaled by the
        incident flux I0*mu0) and ``{"reflectance","transmittance"}`` (the
        normalised ratios R and alpha). Keyed as ``data[mu0][quantity]``.
    """
    tau = np.asarray(tau, dtype=float)
    if leg_coeffs is None:
        leg_coeffs = henyey_greenstein_legendre(n_layers=max(len(tau), 1))

    out = {}
    for mu0 in mu0_values:
        refl = np.zeros_like(tau)
        tdiff = np.zeros_like(tau)
        tdir = np.zeros_like(tau)
        for i, tau_i in enumerate(tau):
            # only_flux=False matches the original (full intensity solve).
            refl[i], tdiff[i], tdir[i] = _slab_fluxes(
                tau_i, ssa, mu0, leg_coeffs[i, :], only_flux=False)
        scale = I0 * mu0
        out[mu0] = {
            "F_u": np.real(refl * scale),
            "F_d": np.real(tdiff * scale),
            "F_s": np.real(tdir * scale),
            "reflectance": np.real(refl),
            "transmittance": np.real(tdiff + tdir),
        }
    return out


def alloy_transmittance_database(alloy_list, elemental_data, L_um_list,
                                 ssa=0.5, mu0=1.0, n_aly=2.00, leg_coeffs=None):
    """Transmittance spectra for many alloys at several thicknesses (Fig. S6).

    For each alloy the minimum-MFP composition is found (BvK + Tamura), the
    optical-depth spectrum is built for each thickness, restricted to the
    trusted tau range, and solved with DISORT at normal incidence.

    Parameters
    ----------
    alloy_list : sequence of str
        Alloy names "ElementA ElementB".
    elemental_data : dict
        Per-element property table.
    L_um_list : sequence of float
        Scattering-layer thicknesses (micrometres).
    ssa : float
        Single-scattering albedo (default 0.5).
    mu0 : float
        Incidence cosine (default 1.0).
    n_aly : float
        Tamura frequency exponent (2.00 for this survey).
    leg_coeffs : ndarray, optional
        Phase-function coefficients; default Henyey-Greenstein.

    Returns
    -------
    dict
        ``{alloy: {"min_x_a","min_x_b", L_um: {"nu","transmittance","reflectance"}}}``
        with frequency ``nu`` in Hz.
    """
    if leg_coeffs is None:
        leg_coeffs = henyey_greenstein_legendre(n_layers=100)

    out = {}
    for alloy in alloy_list:
        element_a, element_b = alloy.split()
        if element_a == element_b:
            continue
        omega, _tau_life, mfp, min_x_a, min_x_b = scattering.alloy_mfp_min_composition(
            element_a, element_b, elemental_data, n_aly=n_aly)
        out[alloy] = {"min_x_a": min_x_a, "min_x_b": min_x_b}

        for L_um in L_um_list:
            tau_acou = L_um * 1e-6 / mfp
            tau_lim, nu_lim = dataio.filter_acoustic_thickness(tau_acou, omega, 1e-10, 5e4)
            refl = np.zeros_like(tau_lim)
            tdiff = np.zeros_like(tau_lim)
            tdir = np.zeros_like(tau_lim)
            for i in range(len(tau_lim)):
                refl[i], tdiff[i], tdir[i] = _slab_fluxes(
                    tau_lim[i], ssa, mu0, leg_coeffs[i, :])
            out[alloy][L_um] = {
                "nu": nu_lim,
                "transmittance": np.real(tdiff + tdir),
                "reflectance": np.real(refl),
            }
    return out

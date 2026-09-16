#!/usr/bin/env python
"""
Regenerate every manuscript figure from scratch.

Runs the full pipeline:
  DFT DOS  ->  mass-difference MFP  ->  DISORT transport  +  AMM interfaces
  ->  Figs 2, 3, 4, 5, 6, S1, S4, S5, S6, S7, S8.

Outputs are written to ``figures/``. This script mirrors, in plain-script form,
the walkthrough in ``notebooks/reproduce_figures.ipynb``.

Usage:
    python generate_all_figures.py
"""

import os
os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np

import phonon_filters as pf
from phonon_filters import dataio, scattering, transport, figures

OUT = "figures"


def main():
    print("Loading DFT density of states and computing SiGe mean free path ...")
    nu_THz, _, dos = dataio.read_dos_average(dataio.default_dos_paths(), num_points=110)
    _, mfp = scattering.mass_difference_mfp_dft(nu_THz, dos)
    nu_tr_Hz, mfp_tr, _ = dataio.truncate_mfp(nu_THz * 1e12, mfp, dos)
    nu_GHz = nu_tr_Hz / 1e9
    elemental = pf.load_elemental_properties()

    # ---- Fig. 2: cumulative transmittance (both stack directions) ----
    print("Fig 2: cumulative transmittance (Si/SiGe/Al) ...")
    mu0_array = np.linspace(0.01, 1.0, 21)
    scat_matrix = transport.transmittance_matrix(mfp_tr, mu0_array, L_um=10, ssa=1 - 1e-6)
    figures.figure2(["Si", "SiGe_avg", "Al"], nu_GHz, mu0_array, scat_matrix,
                    pf.MATERIALS, shade_high_freq=True, outdir=OUT)
    figures.figure2(["Al", "SiGe_avg", "Si"], nu_GHz, mu0_array, scat_matrix,
                    pf.MATERIALS, elemental_data=elemental, temp_overlay=True, outdir=OUT)

    # ---- Fig. 3: parameter sensitivity ----
    print("Fig 3: transmittance sensitivity ...")
    base_params = {"m_Si": 28.08 * 1.66e-27, "m_Ge": 72.08 * 1.66e-27, "V0": 48e-30,
                   "v_La": 6281, "v_Ta": 2937, "n_aly": 2, "DoS": dos}
    param_factors = {"m_Si": 0.25, "m_Ge": 1.75, "v_La": 1.75, "v_Ta": 1.75, "DoS": 1.75}
    lengths = [0.1, 1.0, 10.0]
    mu0_list = [mu0_array[-1], mu0_array[15], mu0_array[6]]
    stack, nu3 = figures.compute_sensitivity(
        nu_THz, dos, pf.MATERIALS, base_params, param_factors,
        lengths, mu0_list, ["Si", "SiGe_avg", "Al"], ssa=0.5)
    figures.figure3(stack, nu3, mu0=1.0, lengths=lengths, param_factors=param_factors, outdir=OUT)

    # ---- Figs. 4 / S8: albedo sweeps ----
    print("Figs 4 & S8: single-scattering-albedo sweep ...")
    albedos = [0.001, 0.5, 0.99, 0.999999]
    sweep = transport.flux_albedo_sweep(mfp_tr, albedos, mu0=1.0, L_um=1.0)
    figures.figure4(nu_GHz, sweep, albedos, outdir=OUT)
    figures.figureS8(nu_GHz, sweep, albedos, outdir=OUT)

    # ---- Figs. 5 / 6 / S1: single-interface AMM ----
    print("Figs 5, 6, S1: single-interface acoustic mismatch ...")
    figures.figure5(pf.MATERIALS, outdir=OUT)
    figures.figure6(pf.MATERIALS, outdir=OUT)
    figures.figureS1(pf.MATERIALS, outdir=OUT)

    # ---- Fig. S4: alloy MFP vs composition ----
    print("Fig S4: alloy MFP vs composition ...")
    figures.figureS4("Silicon", "Germanium", elemental, outdir=OUT)

    # ---- Fig. S5: Debye integrand ----
    print("Fig S5: Debye specific-heat integrand ...")
    figures.figureS5(x_max=25.0, outdir=OUT)

    # ---- Fig. S6: broad alloy survey ----
    print("Fig S6: alloy transmittance survey (this is the slow step) ...")
    database = transport.alloy_transmittance_database(
        pf.FULLY_SOLUBLE_ALLOYS, elemental, [0.1, 1.0, 10.0], ssa=0.5, n_aly=2.00)
    figures.figureS6(database, pf.FULLY_SOLUBLE_ALLOYS, [0.1, 1.0, 10.0], outdir=OUT)

    # ---- Fig. S7: component fluxes vs optical depth ----
    print("Fig S7: component fluxes vs optical depth ...")
    tau = np.logspace(-3, 3, 100)
    fluxes = transport.flux_vs_thickness(tau, (1.0, 0.1), ssa=0.999999)
    figures.figureS7(tau, fluxes, outdir=OUT)

    print(f"Done. Figures written to {OUT}/")


if __name__ == "__main__":
    main()

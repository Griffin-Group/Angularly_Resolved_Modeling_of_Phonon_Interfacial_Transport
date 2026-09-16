"""
Manuscript figure generators.

One function per figure in "Theoretical Design of Phonon Filters for Quantum
Sensing and Qubits". Each function computes the figure data from the physics
modules and renders the plot, returning the underlying arrays so results can be
inspected or reused. Figures are written to ``outdir`` (default ``figures/``).

Figure map
----------
    figure2  : Fig. 2a/2b  cumulative transmittance, Si/SiGe/Al (both directions)
    figure3  : Fig. 3      transmittance sensitivity to model parameters
    figure4  : Fig. 4      transmittance vs frequency for several albedos
    figure5  : Fig. 5      single-interface transmission vs angle (4 panels)
    figure6  : Fig. 6      single-interface transmission vs angle, Nb (2 panels)
    figureS4 : Fig. S5     alloy MFP vs composition
    figureS5 : Fig. S6     Debye specific-heat integrand + cumulative integral
    figureS6 : Fig. S7     transmittance vs frequency for many alloys
    figureS7 : Fig. S8     component fluxes vs optical depth
    figureS8 : Fig. S9     reflectance vs frequency for several albedos
"""

from __future__ import annotations

import copy
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc_context
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D

from . import amm, scattering, thermal, transport

_OUTDIR = Path(__file__).resolve().parent.parent / "figures"
_RC = {"figure.dpi": 300, "font.size": 18, "savefig.bbox": "tight", "font.family": "Arial"}

# Default figure output format. PDF (vector) by default -- publication-ready and
# resolution-independent. Raster PNG output is one uncomment away in ``_savefig``
# (kept, not deleted, so flipping back is trivial). Set _FIG_EXT = "png" to make
# raster the default again. (Sean directive 2026-06-05.)
_FIG_EXT = "png"

# Temperatures (K) for the thermal-occupation overlay in Fig. 2b.
_OVERLAY_TEMPS = (0.05, 0.5, 5.0)

# Aluminium pair-breaking gap frequency 2*Delta/h (GHz). Onset of the
# harmful-quasiparticle band shaded in Fig. 2a; matches the "Al" reference line
# in Fig. S6 and the manuscript Fig. 2a caption ("nu_p >= 2*Delta/h for Al").
_AL_GAP_GHZ = 85.0


def _outpath(outdir, name):
    outdir = Path(_OUTDIR if outdir is None else outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    return outdir / name


def _savefig(fig, outdir, stem):
    """Save ``fig`` as ``<stem>.<_FIG_EXT>`` in ``outdir`` (PDF by default).

    Centralizes the output-format choice in one place so the whole repo flips
    PDF<->PNG by editing ``_FIG_EXT`` (or uncommenting the raster line below)
    rather than touching every figure generator. The PNG ``savefig`` line is
    preserved-but-commented per Sean's directive: uncomment it to ALSO emit a
    300-dpi raster alongside the vector PDF.
    """
    fig.savefig(_outpath(outdir, f"{stem}.{_FIG_EXT}"))
    # fig.savefig(_outpath(outdir, f"{stem}.png"), dpi=300)  # raster PNG (uncomment to also emit)


# ======================================================================
# Fig. 2 -- cumulative transmittance (AMM x mass-difference scattering)
# ======================================================================
def figure2(stack_keys, nu_lim_GHz, mu0_array, scat_matrix, materials,
            elemental_data=None, temp_overlay=False, shade_high_freq=False,
            outdir=None, show=False, show_colorbar=True):

    angles_rad = np.arccos(mu0_array)
    freq = np.asarray(nu_lim_GHz)
    stack = [materials[k] for k in stack_keys]

    # Original behavior (NO reordering fix applied)
    scat = scat_matrix[::-1, :]

    results = {}
    for mode in ("L", "T"):
        t_grid, _ = amm.stack_transmission_isotropic(
            angles_rad, stack, mode=mode
        )
        results[mode] = np.tile(t_grid[:, None], (1, len(freq))) * scat

    mu = np.cos(np.sort(angles_rad))

    with rc_context(rc=_RC):
        fig, ax = plt.subplots(figsize=(10, 8))
        cmap = plt.get_cmap("rainbow")
        idx = np.linspace(0, len(mu) - 1, 15, dtype=int)
        norm = mcolors.Normalize(vmin=mu[idx].min(), vmax=mu[idx].max())

        if temp_overlay:
            ax2 = ax.twinx()
            ax2.set_ylabel(r"Normalized $DoS \cdot f_{BE}$")
            from .dataio import element_symbol_to_name
            src_name = element_symbol_to_name(stack_keys[0]) or "Aluminium"
            alphas = np.linspace(1.0, 0.3, len(_OVERLAY_TEMPS))
            for a, T in zip(alphas, _OVERLAY_TEMPS):
                _, _, dos_n, _ = thermal.dos_times_occupation(
                    src_name,
                    nu_lim_GHz * 2 * np.pi * 1e9,
                    T,
                    elemental_data
                )
                ax2.semilogx(
                    nu_lim_GHz,
                    dos_n / np.max(dos_n),
                    color="black",
                    alpha=a,
                    label=f"{T} K"
                )

        for j in idx:
            color = cmap(norm(mu[j]))
            ax.semilogx(freq, results["L"][j, :], color=color, linestyle="-",
                        label="LA" if j == idx[0] else None)
            ax.semilogx(freq, results["T"][j, :], color=color, linestyle=":",
                        label="TA" if j == idx[0] else None)

        if shade_high_freq:
            ax.axvspan(_AL_GAP_GHZ, freq.max(), color="blue", alpha=0.075)

        if show_colorbar:
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            fig.colorbar(sm, ax=ax).set_label(r"$\mu = \cos(\theta)$")

        ax.set_xlabel(r"Phonon frequency, $\nu_p$ (GHz)")
        ax.set_ylabel(r"Transmittance, $\alpha_{total}$")
        ax.grid(True)

        if temp_overlay:
            h1, l1 = ax.get_legend_handles_labels()
            h2, l2 = ax2.get_legend_handles_labels()
            ax.legend(h1 + h2, l1 + l2, loc="best")
        else:
            ax.legend(loc="upper right")

        fig.tight_layout()
        _savefig(fig, outdir, "fig2_" + "_".join(stack_keys))

        plt.show() if show else plt.close(fig)

    return {
        "L": results["L"],
        "T": results["T"],
        "mu": mu,
        "freq_GHz": freq
    }


# ======================================================================
# Fig. 3 -- parameter sensitivity of the cumulative transmittance
# ======================================================================
def compute_sensitivity(nu_THz, dos, materials, base_params, param_factors,
                        lengths, mu0_list, stack_keys, Mode, ssa=0.5, 
                        consistent_velocity_scaling=True, ):

    stack_out = {case: {} for case in ["base"] + list(param_factors)}
    nu_trunc_GHz = None

    mu0_list = np.sort(mu0_list)

    for L_um in lengths:
        L_m = L_um * 1e-6

        for case in ["base"] + list(param_factors):
            params = copy.deepcopy(base_params)
            mats = copy.deepcopy(materials)

            if case != "base":
                params[case] *= param_factors[case]

            if consistent_velocity_scaling:
                if case == "v_La":
                    mats["SiGe_avg"]["cL"] *= param_factors[case]
                elif case == "v_Ta":
                    mats["SiGe_avg"]["cT"] *= param_factors[case]

            _, mfp = scattering.mass_difference_mfp_dft(
                nu_THz, params["DoS"],
                m_Si=params["m_Si"],
                m_Ge=params["m_Ge"],
                V0=params["V0"],
                v_La=params["v_La"],
                v_Ta=params["v_Ta"]
            )

            mfp = np.clip(mfp, 1e-30, None)

            from .dataio import truncate_mfp
            nu_tr, mfp_tr, _ = truncate_mfp(
                nu_THz * 1e12, mfp, params["DoS"]
            )
            nu_trunc_GHz = nu_tr / 1e9

            tau_scaled = L_m / mfp_tr

            leg = transport.henyey_greenstein_legendre(
                n_layers=max(len(tau_scaled), 1)
            )

            t_all = []

            for mu0 in mu0_list:
                tdiff = np.zeros_like(tau_scaled)
                tdir = np.zeros_like(tau_scaled)

                for i in range(len(tau_scaled)):
                    _, td, tdr = transport._slab_fluxes(
                        tau_scaled[i],
                        ssa,
                        mu0,
                        leg[i, :]
                    )
                    tdiff[i], tdir[i] = td, tdr

                t_all.append(np.real(tdiff + tdir))

            t_all = np.vstack(t_all)
            t_all = t_all[::-1]

            mode = Mode
            stack = [mats[k] for k in stack_keys]


            t_grid, _ = amm.stack_transmission_isotropic(
                np.arccos(mu0_list),
                stack,
                mode=mode
            )

            alpha = t_grid[:, None] * t_all

            stack_out[case][L_um] = {
                mu0: alpha[k, :]
                for k, mu0 in enumerate(mu0_list[::-1])
            }

    return stack_out, nu_trunc_GHz


def figure3(stack, nu_trunc_GHz, mu0, lengths, param_factors,
            highlight_length=1.0, outdir=None, show=False):
    """Plot the transmittance sensitivity (Fig. 3).

    Parameters
    ----------
    stack : dict
        Output of :func:`compute_sensitivity`.
    nu_trunc_GHz : ndarray
        Frequency grid (GHz).
    mu0 : float
        Incidence cosine to plot (e.g. 1.0).
    lengths : sequence of float
        Thicknesses (micrometres); each gets a distinct line style.
    param_factors : dict
        Parameter -> scale factor (defines colours and legend).
    highlight_length : float
        Thickness drawn with the thick base curve.
    """
    linestyles = ["--", "-", ":"]
    colors = plt.get_cmap("rainbow")(np.linspace(0, 1, len(param_factors)))
    with rc_context(rc=_RC):
        fig, ax = plt.subplots(figsize=(9, 7))
        for j, L_um in enumerate(lengths):
            for p, param in enumerate(param_factors):
                ax.semilogx(nu_trunc_GHz, np.real(stack[param][L_um][mu0]),
                            color=colors[p], linestyle=linestyles[j], lw=2, alpha=0.5)
            lw = 3.5 if L_um == highlight_length else 2
            ax.semilogx(nu_trunc_GHz, np.real(stack["base"][L_um][mu0]),
                        color="black", linestyle=linestyles[j], lw=lw)

        legend = [Line2D([0], [0], color="black", lw=3, label="Base case")]
        for p, param in enumerate(param_factors):
            if "_" in param:
                b, s = param.split("_", 1)
                lbl = rf"${b}_{{{s}}}$"
            else:
                lbl = rf"${param}$"
            legend.append(Line2D([0], [0], color=colors[p], lw=3, alpha=0.5,
                                 label=f"{lbl} (x{param_factors[param]:.2f})"))
        for j, L_um in enumerate(lengths):
            legend.append(Line2D([0], [0], color="black", linestyle=linestyles[j],
                                 lw=3, label=f"L={L_um} μm"))
        ax.legend(handles=legend, fontsize=12)
        ax.set_xlabel(r"Phonon frequency, $\nu_p$ (GHz)")
        ax.set_ylabel(r"Transmittance, $\alpha_{total}$")
        ax.grid(True, which="both", linestyle="--", linewidth=0.5)
        ax.set_xlim(50, 2e3)
        fig.tight_layout()
        _savefig(fig, outdir, "fig3_sensitivity")
        plt.show() if show else plt.close(fig)


# ======================================================================
# Fig. 4 / S8 -- albedo sweeps
# ======================================================================
_ALBEDO_COLORS = {0.001: "red", 0.5: "black", 0.99: "blue", 0.999999: "green"}


def _albedo_plot(nu_GHz, sweep, albedos, quantity, ylabel, stem, outdir, show):
    with rc_context(rc=_RC):
        fig, ax = plt.subplots(figsize=(9, 7))
        legend = []
        for ssa in albedos:
            color = _ALBEDO_COLORS.get(ssa, "gray")
            lw = 3 if ssa == 0.5 else 2
            ax.semilogx(nu_GHz, sweep[ssa][quantity], color=color, linewidth=lw)
            legend.append(Line2D([0], [0], color=color, lw=3, label=rf"$\sigma_p$ = {ssa}"))
        ax.legend(handles=legend, title="Single Scattering Albedo")
        ax.set_xlabel(r"Phonon frequency, $\nu_p$ (GHz)")
        ax.set_ylabel(ylabel)
        ax.grid(True, which="both", linestyle="--", linewidth=0.5)
        fig.tight_layout()
        ax.set_xlim(50, 2e3)
        _savefig(fig, outdir, stem)
        plt.show() if show else plt.close(fig)


def figure4(nu_trunc_GHz, sweep, albedos, outdir=None, show=False):
    """Transmittance vs frequency for several single-scattering albedos (Fig. 4).

    Parameters
    ----------
    nu_trunc_GHz : ndarray
        Frequency grid (GHz).
    sweep : dict
        Output of ``transport.flux_albedo_sweep``.
    albedos : sequence of float
        Albedos to plot, in legend order.
    """
    _albedo_plot(nu_trunc_GHz, sweep, albedos, "transmittance",
                 r"Transmittance, $\alpha_{mass-difference}$",
                 "fig4_transmittance_albedo", outdir, show)


def figureS9(nu_trunc_GHz, sweep, albedos, outdir=None, show=False):
    """Reflectance vs frequency for several single-scattering albedos (Fig. S8)."""
    _albedo_plot(nu_trunc_GHz, sweep, albedos, "reflectance",
                 r"Reflectance, $R_{mass-difference}$",
                 "figS8_reflectance_albedo", outdir, show)


# ======================================================================
# Fig. 5 / 6 -- single-interface angular transmission
# ======================================================================
def _interface_panel(ax, mat1, mat2, materials, title):
    mu = np.linspace(0.0, 1.0, 500)
    aL = amm.interface_power_transmission(materials[mat1], materials[mat2], "L", mu)
    aT = amm.interface_power_transmission(materials[mat1], materials[mat2], "T", mu)
    ax.plot(mu, aT, color="tab:blue", linestyle="-", label="TA")
    ax.plot(mu, aL, color="tab:orange", linestyle="--", label="LA")
    ax.set_xlabel(r"$\cos\,\Theta$")
    ax.set_ylabel(r"$\alpha$")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.05)
    ax.set_title(title, fontsize=12)
    ax.legend(loc="lower right", fontsize=10)
    ax.grid(True, alpha=0.3)


def figure5(materials, outdir=None, show=False):
    """Single-interface transmission vs angle for the Al-qubit stack (Fig. 5).

    Four interfaces (Si->SiGe, SiGe->Al, Al->SiGe, SiGe->Si), TA solid and LA
    dashed. The critical-angle cutoffs follow from Snell refraction into the
    faster medium.
    """
    xticks = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    panels = [("SiGe_avg", "Al", "(a) SiGe → Al"),
              ("SiGe_avg", "Nb", "(b) SiGe → Nb"),
              ("SiGe_avg", "NbN", "(c) SiGe → NbN"),
              ("Al", "SiGe_avg", "(d) Al → SiGe"),
              ("Nb", "SiGe_avg", "(e) Nb → SiGe"),
              ("NbN", "SiGe_avg", "(f) NbN → SiGe"),]
    with rc_context(rc=_RC):
        fig, axes = plt.subplots(2, 3, figsize=(15, 9))
        for ax, (m1, m2, title) in zip(axes.flat, panels):
            _interface_panel(ax, m1, m2, materials, title)
            ax.set_xticks(xticks)
        fig.tight_layout()
        _savefig(fig, outdir, "fig5_interface_transmission")
        plt.show() if show else plt.close(fig)


# ======================================================================
# Fig. S5 -- alloy MFP vs composition
# ======================================================================
def figureS5(element_a, element_b, elemental_data, freqs_THz=(0.1, 0.7, 1.1),
             n_aly=2.047, outdir=None, show=False):
    """Composition-dependent alloy mean free path (Fig. S4)."""
    from .dataio import element_name_to_symbol
    x_b, mfp_um, min_x_a, min_x_b = scattering.mfp_vs_composition(
        element_a, element_b, elemental_data, freqs_THz=freqs_THz, n_aly=n_aly)
    sym_a = element_name_to_symbol(element_a)
    sym_b = element_name_to_symbol(element_b)
    colors = ["tab:blue", "tab:orange", "tab:green"]
    with rc_context(rc={**_RC, "font.size": 14}):
        fig, ax = plt.subplots(figsize=(8, 6))
        for i, f_THz in enumerate(freqs_THz):
            ax.loglog(x_b, mfp_um[i], color=colors[i % len(colors)],
                      label=rf"$\nu_p = {f_THz * 1e3:.0f}\,\mathrm{{GHz}}$")
        ax.set_title(rf"MFP vs concentration | $n_{{aly}} = {n_aly}$")
        ax.set_xlabel(rf"$x$ in ${sym_a}_{{1-x}}{sym_b}_x$")
        ax.set_ylabel(r"Mean free path ($\mu$m)")
        ax.legend()
        ax.grid(True, which="both", linestyle="--", linewidth=0.5)
        fig.tight_layout()
        _savefig(fig, outdir, "figS4_mfp_composition")
        plt.show() if show else plt.close(fig)
    return {"x_b": x_b, "mfp_um": mfp_um, "min_x_a": min_x_a, "min_x_b": min_x_b}


# ======================================================================
# Fig. S6 -- Debye specific-heat integrand
# ======================================================================
def figureS6(x_max=15.0, outdir=None, show=False):
    """Debye specific-heat integrand and its cumulative integral (Fig. S5)."""
    x, integrand, cumulative = thermal.debye_integrand_cumulative(x_max=x_max)
    with rc_context(rc={**_RC, "font.size": 14}):
        fig, ax1 = plt.subplots(figsize=(7, 5))
        ax1.plot(x, integrand, color="blue", lw=2)
        ax1.set_xlabel(r"$x = h\nu / (k_B T_q)$")
        ax1.set_ylabel(r"$\frac{x^4 e^x}{(e^x-1)^2}$", color="blue")
        ax1.tick_params(axis="y", labelcolor="blue")
        ax2 = ax1.twinx()
        ax2.plot(x, cumulative, color="red", lw=2, ls="--")
        ax2.set_ylabel(r"$\int_0^x \frac{t^4 e^t}{(e^t-1)^2}\,dt$", color="red")
        ax2.tick_params(axis="y", labelcolor="red")
        ax1.grid(alpha=0.3)
        fig.tight_layout()
        _savefig(fig, outdir, "figS5_debye_integrand")
        plt.show() if show else plt.close(fig)
    return {"x": x, "integrand": integrand, "cumulative": cumulative}


# ======================================================================
# Fig. S7 -- transmittance vs frequency for many alloys
# ======================================================================
def figureS7(database, alloys, L_um_list, group_size=4, outdir=None, show=False):
    """Transmittance spectra for a survey of solid-solution alloys (Fig. S6).

    Parameters
    ----------
    database : dict
        Output of ``transport.alloy_transmittance_database``.
    alloys : sequence of str
        Alloy names to plot (subset of database keys).
    L_um_list : sequence of float
        Thicknesses; each gets a distinct line style.
    group_size : int
        Alloys per subplot.
    """
    from .dataio import element_name_to_symbol
    alloys = [a for a in alloys if a in database]
    groups = [alloys[i:i + group_size] for i in range(0, len(alloys), group_size)]
    nrows = int(np.ceil(len(groups) / 2))
    linestyles = ["--", "-", ":", "-."]
    colors = plt.cm.tab10.colors

    with rc_context(rc={**_RC, "font.size": 14}):
        fig, axes = plt.subplots(nrows, 2, figsize=(10, 4 * nrows), sharey=True)
        axes = np.atleast_1d(axes).flatten()
        for gi, group in enumerate(groups):
            ax = axes[gi]
            for ci, alloy in enumerate(group):
                a, b = alloy.split()
                sym = (element_name_to_symbol(a) or a) + (element_name_to_symbol(b) or b)
                color = colors[ci % len(colors)]
                for li, L in enumerate(L_um_list):
                    if L not in database[alloy]:
                        continue
                    d = database[alloy][L]
                    ax.semilogx(d["nu"] / 1e9, d["transmittance"], color=color,
                                linestyle=linestyles[li % len(linestyles)], lw=2)
                ax.plot([], [], color=color, label=sym)
            # Characteristic-frequency reference lines: the TA thermal cutoff at
            # 50 mK and the pair-breaking gap frequencies 2*Delta/h for Al, Nb,
            # NbN. The text labels are placed vertically at the top of each line
            # (data x, axes-fraction y via get_xaxis_transform) -- previously the
            # `txt` was computed but never rendered.
            for nu_ref, txt in [(85, "Al"),
                                (750, "Nb"), (1693, "NbN")]:
                ax.axvline(nu_ref, color="gray", linestyle="--", linewidth=1)
                ax.text(nu_ref, 0.5, txt, transform=ax.get_xaxis_transform(),
                        rotation=90, va="top", ha="right", fontsize=8, color="black")
            ax.set_xlabel(r"$\nu_p$ (GHz)")
            ax.set_ylabel(r"Transmittance, $\alpha_{mass-diff}$")
            ax.grid(True, which="both", ls="--", alpha=0.5)
            ax.set_xlim(50,None)
            ax.legend(fontsize=9, loc="best")
            
        for j in range(len(groups), len(axes)):
            axes[j].axis("off")
        length_handles = [Line2D([0], [0], color="black", linestyle=ls, lw=2)
                          for ls in linestyles[:len(L_um_list)]]
        fig.legend(length_handles, [rf"$L={L}\,\mu$m" for L in L_um_list],
                   loc="center", bbox_to_anchor = (0.5, -0.025), ncol=len(L_um_list), frameon=False)
        fig.tight_layout()
        _savefig(fig, outdir, "figS6_alloy_survey")
        plt.show() if show else plt.close(fig)


# ======================================================================
# Fig. S8 -- component fluxes vs optical depth
# ======================================================================
_FLUX_STYLE = {
    "F_u": ("orange", r"$F_u$"),
    "F_d": ("red", r"$F_d$"),
    "F_s": ("violet", r"$F_s$"),
    "reflectance": ("cyan", r"$R=F_u/F_i$"),
    "transmittance": ("blue", r"$\alpha=(F_s+F_d)/F_i$"),
}


def figureS8(tau, fluxes, quantities=("F_u", "F_d", "F_s", "reflectance", "transmittance"),
             ssa=0.999999, g=transport.ASYMMETRY_G, outdir=None, show=False):
    """Component fluxes versus optical depth (Fig. S7).

    Parameters
    ----------
    tau : ndarray
        Optical depth grid.
    fluxes : dict
        Output of ``transport.flux_vs_thickness`` (keyed by mu0).
    quantities : sequence of str
        Which curves to draw.
    ssa, g : float
        Albedo and HG asymmetry, for the title.
    """
    with rc_context(rc={**_RC, "font.size": 14}):
        fig, ax = plt.subplots(figsize=(7, 6))
        mu0_values = list(fluxes)
        for i, mu0 in enumerate(mu0_values):
            ls = "-" if i == 0 else "--"
            for q in quantities:
                color, _ = _FLUX_STYLE[q]
                ax.loglog(tau, fluxes[mu0][q], linestyle=ls, color=color, lw=2)
        handles = [Line2D([0], [0], color=_FLUX_STYLE[q][0], label=_FLUX_STYLE[q][1])
                   for q in quantities]
        for i, mu0 in enumerate(mu0_values):
            handles.append(Line2D([0], [0], color="green",
                                  linestyle="-" if i == 0 else "--", label=rf"$\mu_0={mu0}$"))
        ax.legend(handles=handles, loc="upper left", fontsize=10)
        ax.set_ylim(1e-3, 10)
        ax.set_xlabel(r"$\tau_p$")
        ax.set_ylabel(r"Component fluxes ($F_s, F_u, F_d$)")
        ax.set_title("Different fluxes vs acoustic thickness\n"
                     rf"($F_i=\pi$, g={g}, $\sigma_p$={ssa})")
        ax.grid()
        fig.tight_layout()
        _savefig(fig, outdir, "figS7_fluxes")
        plt.show() if show else plt.close(fig)

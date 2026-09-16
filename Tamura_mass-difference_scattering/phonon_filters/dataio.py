"""
Input/output and array-conditioning helpers.

Covers:
  * reading the DFT phonon density-of-states (DOS) files and averaging across
    SiGe phases onto a common frequency grid,
  * converting between element names and chemical symbols (for plot labels),
  * truncating / filtering frequency-dependent arrays to a physical range.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
from scipy.interpolate import interp1d

# Canonical element name -> chemical symbol (used only for figure labels).
NAME_TO_SYMBOL: dict[str, str] = {
    "Hydrogen": "H", "Helium": "He", "Lithium": "Li", "Beryllium": "Be", "Boron": "B", "Carbon": "C",
    "Nitrogen": "N", "Oxygen": "O", "Fluorine": "F", "Neon": "Ne", "Sodium": "Na", "Magnesium": "Mg",
    "Aluminium": "Al", "Silicon": "Si", "Phosphorus": "P", "Sulfur": "S", "Chlorine": "Cl", "Argon": "Ar",
    "Potassium": "K", "Calcium": "Ca", "Scandium": "Sc", "Titanium": "Ti", "Vanadium": "V", "Chromium": "Cr",
    "Manganese": "Mn", "Iron": "Fe", "Cobalt": "Co", "Nickel": "Ni", "Copper": "Cu", "Zinc": "Zn",
    "Gallium": "Ga", "Germanium": "Ge", "Arsenic": "As", "Selenium": "Se", "Bromine": "Br", "Krypton": "Kr",
    "Rubidium": "Rb", "Strontium": "Sr", "Yttrium": "Y", "Zirconium": "Zr", "Niobium": "Nb", "Molybdenum": "Mo",
    "Technetium": "Tc", "Ruthenium": "Ru", "Rhodium": "Rh", "Palladium": "Pd", "Silver": "Ag", "Cadmium": "Cd",
    "Indium": "In", "Tin": "Sn", "Antimony": "Sb", "Tellurium": "Te", "Iodine": "I", "Xenon": "Xe",
    "Cesium": "Cs", "Barium": "Ba", "Lanthanum": "La", "Cerium": "Ce", "Praseodymium": "Pr", "Neodymium": "Nd",
    "Promethium": "Pm", "Samarium": "Sm", "Europium": "Eu", "Gadolinium": "Gd", "Terbium": "Tb", "Dysprosium": "Dy",
    "Holmium": "Ho", "Erbium": "Er", "Thulium": "Tm", "Ytterbium": "Yb", "Lutetium": "Lu", "Hafnium": "Hf",
    "Tantalum": "Ta", "Tungsten": "W", "Rhenium": "Re", "Osmium": "Os", "Iridium": "Ir", "Platinum": "Pt",
    "Gold": "Au", "Mercury": "Hg", "Thallium": "Tl", "Lead": "Pb", "Bismuth": "Bi", "Polonium": "Po",
    "Astatine": "At", "Radon": "Rn", "Francium": "Fr", "Radium": "Ra", "Actinium": "Ac", "Thorium": "Th",
    "Protactinium": "Pa", "Uranium": "U", "Neptunium": "Np", "Plutonium": "Pu", "Americium": "Am", "Curium": "Cm",
    "Berkelium": "Bk", "Californium": "Cf", "Einsteinium": "Es", "Fermium": "Fm", "Mendelevium": "Md", "Nobelium": "No",
    "Lawrencium": "Lr", "Rutherfordium": "Rf", "Dubnium": "Db", "Seaborgium": "Sg", "Bohrium": "Bh", "Hassium": "Hs",
    "Meitnerium": "Mt", "Darmstadtium": "Ds", "Roentgenium": "Rg", "Copernicium": "Cn", "Nihonium": "Nh", "Flerovium": "Fl",
    "Moscovium": "Mc", "Livermorium": "Lv", "Tennessine": "Ts", "Oganesson": "Og",
}
SYMBOL_TO_NAME: dict[str, str] = {sym: name for name, sym in NAME_TO_SYMBOL.items()}


def element_name_to_symbol(name: str) -> str | None:
    """'Silicon' -> 'Si'. Returns None if unknown."""
    return NAME_TO_SYMBOL.get(name.strip().capitalize())


def element_symbol_to_name(symbol: str) -> str | None:
    """'Si' -> 'Silicon'. Returns None if unknown."""
    return SYMBOL_TO_NAME.get(symbol.strip().capitalize())


def round_sig(x: float, sig: int) -> float:
    """Round ``x`` to ``sig`` significant figures."""
    if x == 0:
        return 0.0
    return round(x, sig - int(math.floor(math.log10(abs(x)))) - 1)


def read_dos_average(filenames, x_start: float = 0.005, num_points: int = 110):
    """Read several two-column DOS files and average them on a common grid.

    Each file is whitespace-delimited ``frequency  DOS`` with optional ``#``
    comment lines. Negative or zero frequencies (acoustic-sum-rule noise near
    Gamma) are dropped. Every dataset is linearly interpolated onto one shared
    log-spaced frequency grid spanning ``x_start`` to the largest frequency
    found across all files, and the per-grid-point mean DOS is returned. This
    produces the phase-averaged SiGe DOS used in the main text.

    Parameters
    ----------
    filenames : sequence of str or Path
        Paths to the DOS files (frequency in THz, DOS per unit cell).
    x_start : float
        Lower end of the common frequency grid (THz). Default 0.005.
    num_points : int
        Number of points on the common log-spaced grid. Default 110.

    Returns
    -------
    freq_common : ndarray
        Common frequency grid (THz).
    dos_all : ndarray, shape (n_files, num_points)
        Interpolated DOS for each input file.
    dos_avg : ndarray
        Mean DOS across files (per unit cell).
    """
    datasets = []
    global_x_max = 0.0
    for filename in filenames:
        x_vals, y_vals = [], []
        with open(filename, "r") as fh:
            for line in fh:
                s = line.strip()
                if not s or s.startswith("#"):
                    continue
                parts = s.split()
                if len(parts) < 2:
                    continue
                try:
                    xv, yv = float(parts[0]), float(parts[1])
                except ValueError:
                    continue
                if xv > 0:
                    x_vals.append(xv)
                    y_vals.append(yv)
        x = np.asarray(x_vals)
        y = np.asarray(y_vals)
        order = np.argsort(x)
        x, y = x[order], y[order]
        datasets.append((x, y))
        global_x_max = max(global_x_max, x.max())

    freq_common = np.logspace(np.log10(x_start), np.log10(global_x_max), num_points)
    dos_all = np.array([
        interp1d(x, y, kind="linear", fill_value="extrapolate")(freq_common)
        for x, y in datasets
    ])
    dos_avg = dos_all.mean(axis=0)
    return freq_common, dos_all, dos_avg


def truncate_mfp(freq_Hz, mfp_m, dos, max_freq: float = 3e12):
    """Keep only finite, non-negative MFP points below ``max_freq``.

    Parameters
    ----------
    freq_Hz, mfp_m, dos : ndarray
        Frequency (Hz), mean free path (m), and DOS arrays of equal length.
    max_freq : float
        Upper frequency cutoff in Hz. Default 3e12 (3 THz).

    Returns
    -------
    freq_Hz, mfp_m, dos : ndarray
        The three arrays restricted to the valid mask.
    """
    mask = (freq_Hz <= max_freq) & np.isfinite(mfp_m) & (mfp_m >= 0)
    return freq_Hz[mask], mfp_m[mask], dos[mask]


def filter_acoustic_thickness(tau, omega, low, high, max_entries: int = 100):
    """Restrict optical depth ``tau`` to ``[low, high]`` and convert omega->Hz.

    Used when mapping a frequency-dependent optical depth onto the DISORT
    solver: values outside the trusted range are dropped and the result is
    down-sampled to at most ``max_entries`` evenly spaced points.

    Parameters
    ----------
    tau : ndarray
        Optical depth (acoustic thickness) per frequency.
    omega : ndarray
        Angular frequency (rad/s) per frequency.
    low, high : float
        Inclusive bounds on ``tau``.
    max_entries : int
        Maximum number of returned points.

    Returns
    -------
    tau_out : ndarray
        Filtered (and possibly down-sampled) optical depth.
    nu_out : ndarray
        Corresponding linear frequency (Hz).
    """
    tau = np.asarray(tau)
    omega = np.asarray(omega)
    mask = (tau >= low) & (tau <= high)
    tau_f = tau[mask]
    nu_f = omega[mask] / (2 * np.pi)
    if len(tau_f) <= max_entries:
        return tau_f, nu_f
    idx = np.linspace(0, len(tau_f) - 1, max_entries, dtype=int)
    return tau_f[idx], nu_f[idx]


# Default DOS files for the phase-averaged SiGe scattering layer (Fig. 2-4).
DEFAULT_SIGE_DOS_FILES = (
    "SiGe_monoclinic_dos.dat",
    "SiGe_hexagonal_dos.dat",
    "SiGe_cubic_dos.dat",
)


def default_dos_paths(data_dir: str | Path | None = None) -> list[Path]:
    """Absolute paths to the three bundled SiGe DOS files."""
    if data_dir is None:
        data_dir = Path(__file__).resolve().parent.parent / "data"
    data_dir = Path(data_dir)
    return [data_dir / name for name in DEFAULT_SIGE_DOS_FILES]

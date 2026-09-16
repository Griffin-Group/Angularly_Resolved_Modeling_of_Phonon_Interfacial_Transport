"""
Material property data for the phonon-filter calculations.

Two kinds of property sets are used by the manuscript:

1. Bulk **acoustic** properties (mass density, longitudinal and transverse
   sound speeds). These feed the Acoustic Mismatch Model (AMM) at the
   interfaces between layers (paper Eq. 2). They live in ``MATERIALS``.

2. Per-element properties (atomic mass, density, polarization-averaged sound
   velocity, atoms per primitive cell). These feed the Born--von Karman (BvK)
   dispersion and the Tamura mass-difference scattering model used for the
   alloy mean free paths (paper Eqs. 4-6). They live in a JSON data file and
   are returned by :func:`load_elemental_properties`.

The ``FULLY_SOLUBLE_ALLOYS`` list is the set of binary alloys (forming solid
solutions over a wide composition range) used for the broad survey in Fig. S6.
"""

from __future__ import annotations

import json
from pathlib import Path

# Acoustic properties used by the AMM interface model.
#   rho : mass density            [kg / m^3]
#   cL  : longitudinal sound speed [m / s]
#   cT  : transverse  sound speed  [m / s]
#
# Two SiGe entries are used by different figures:
#   "SiGe_avg" -- polarization/phase-averaged Si-Ge solid solution; the
#       scattering layer for the cumulative-transmittance results (Figs. 2-4).
#   "SiGe"     -- the F-43m phase values; the SiGe layer used in the angular
#       single-interface AMM results (Figs. 5, 6, S1).
#
# The Nb values below ({4300, 3500}) are those that reproduce the published
# Nb-interface results (Fig. 6, Fig. S1: Nb->SiGe T_eff = 0.866 TA / 0.386 LA).
# A later edit in the original code reverted Nb's cT to 2281; that revert
# postdates the published figures and is not used here (see CORRECTIONS.md).
MATERIALS: dict[str, dict[str, float]] = {
    "Nb":       {"rho": 8527,   "cL": 5080,   "cT": 2340},
    "NbN":      {"rho": 8210,   "cL": 8090,   "cT": 4240}, 
    "Al":       {"rho": 2710.0, "cL": 6860.0, "cT": 3460.0},
    "Si":       {"rho": 2281.0, "cL": 8210.0, "cT": 5730.0},
    "SiGe_avg": {"rho": 3814.0, "cL": 6440.0, "cT": 4160.0},
}

# Binary solid-solution alloys surveyed in Fig. S6. Each entry is
# "ElementA ElementB" with element names matching the keys returned by
# :func:`load_elemental_properties`.
FULLY_SOLUBLE_ALLOYS: list[str] = [
    "Silver Gold",
    "Silver Palladium",
    "Bismuth Antimony",
    "Calcium Strontium",
    "Cobalt Rhenium",
    "Vanadium Chromium",
    "Potassium Rubidium",
    "Molybdenum Niobium",
    "Molybdenum Tantalum",
    "Molybdenum Vanadium",
    "Molybdenum Tungsten",
    "Niobium Tantalum",
    "Niobium Vanadium",
    "Niobium Tungsten",
    "Nickel Platinum",
    "Osmium Rhenium",
    "Osmium Ruthenium",
    "Rhenium Ruthenium",
    "Scandium Yttrium",
    "Selenium Tellurium",
    "Tantalum Tungsten",
    "Titanium Zirconium",
    "Vanadium Tungsten",
    "Silicon Germanium",
]

# Column layout of each value list in elemental_properties.json.
#   [0] atomic_mass            [amu]
#   [1] density                [g / cm^3]
#   [2] sound_velocity         [m / s]   (polarization averaged)
#   [3] atoms_per_primitive_cell
#   [4] melting_point          [K]       (carried for reference; unused here)
ELEMENTAL_COLUMNS = (
    "atomic_mass_amu",
    "density_g_per_cm3",
    "sound_velocity_m_per_s",
    "atoms_per_primitive_cell",
    "melting_point_K",
)

_DATA_DIR = Path(__file__).resolve().parent.parent / "data"


def load_elemental_properties(path: str | Path | None = None) -> dict[str, list]:
    """Return the per-element property table.

    Parameters
    ----------
    path : str or Path, optional
        Location of the JSON file. Defaults to
        ``<repo>/data/elemental_properties.json``.

    Returns
    -------
    dict
        Maps element name -> list of properties in the order described by
        :data:`ELEMENTAL_COLUMNS`.
    """
    if path is None:
        path = _DATA_DIR / "elemental_properties.json"
    with open(path, "r") as fh:
        return json.load(fh)

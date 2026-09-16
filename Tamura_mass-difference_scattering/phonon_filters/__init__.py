"""
phonon_filters
==============

A self-contained reimplementation of the calculations behind

    "Theoretical Design of Phonon Filters for Quantum Sensing and Qubits"
    M. R. Krisshnamurthi, M. A. Hussien, S. Lubner, S. M. Griffin.

The package combines two physical models:

* the **Acoustic Mismatch Model** (AMM) for elastic, specular phonon
  transmission at material interfaces (``amm``), and
* the **Tamura mass-difference scattering** model with **Born--von Karman**
  dispersion for the frequency-dependent mean free path inside the alloy
  scattering layer (``dispersion``, ``scattering``),

solved together through a phonon **radiative-transport** treatment of the
scattering layer (``transport``, via PythonicDISORT). ``figures`` reproduces
every Python-generated manuscript figure.

Module overview
---------------
    materials   : acoustic property tables, alloy lists, elemental data loader
    dataio      : DOS file reading, element<->symbol maps, array truncation
    dispersion  : Born--von Karman group velocity and density of states
    scattering  : Tamura mass-difference lifetimes and mean free paths
    amm         : acoustic-mismatch interface and multilayer transmission
    transport   : DISORT radiative-transport sweeps
    thermal     : Bose-Einstein occupation, DOS weighting, Debye integrand
    figures     : one generator per manuscript figure

See ``notebooks/reproduce_figures.ipynb`` for an end-to-end walkthrough.
"""

from . import (
    materials,
    dataio,
    dispersion,
    scattering,
    amm,
    transport,
    thermal,
    figures,
)

from .materials import MATERIALS, FULLY_SOLUBLE_ALLOYS, load_elemental_properties

__all__ = [
    "materials", "dataio", "dispersion", "scattering", "amm",
    "transport", "thermal", "figures",
    "MATERIALS", "FULLY_SOLUBLE_ALLOYS", "load_elemental_properties",
]

__version__ = "1.0.0"

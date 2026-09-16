# Theoretical Design of Phonon Filters for Quantum Sensing and Qubits

This repository contains the code, input data, calculated results, and
figure-generation materials supporting the manuscript:

> **Theoretical Design of Phonon Filters for Quantum Sensing and Qubits** by Mathan R. K, Musa A. M. Hussien, Sean Lubner, and Sinead M. Griffin

**Manuscript status:** In preparation

The computational workflow computes the frequency-, angle-, and polarization-dependent
phonon transmittance of a Si / SiGe / Al (or Nb or NbN) layer stack, demonstrating how
a SiGe scattering layer acts as a sharp phonon filter that blocks pair-breaking
athermal phonons while passing thermal phonons. Every Python-generated figure in
the manuscript can be regenerated from the bundled inputs with one command.

This repository is a pre-submission research release. The contents may be
updated before or during peer review.


## Physics

Two mechanisms are combined:

1. **Acoustic Mismatch Model (AMM)** — phonons are treated as plane waves that
   refract specularly at interfaces (Snell's law). The interface energy
   transmission is set by the acoustic-impedance contrast,
   `alpha = 1 - R` with `R = ((Z1 cos θ1 - Z2 cos θ2)/(Z1 cos θ1 + Z2 cos θ2))^2`
   (`Z = ρc`). This gives the angular and polarization dependence (Figs. 5, 6, S1).

2. **Tamura mass-difference scattering** — random mass substitution in the SiGe
   alloy scatters phonons at a rate `∝ ω² D(ω)`, yielding a frequency-dependent
   mean free path `Λ(ω)`. The scattering layer is then solved as a phonon
   radiative-transport problem (the phonon analogue of the radiative transfer
   equation) with the discrete-ordinates solver **PythonicDISORT** (D. J. X.
   Ho; `pip install PythonicDISORT`) to obtain the transmittance spectrum
   (Figs. 2, 3, 4, S6, S7, S8).

The cumulative transmittance is the product of the AMM (angular) and
scattering (spectral & angular) contributions.

---
## Installation

### Requirements

- Python 3.11
- NumPy
- Matplotlib
- Phonopy
- scipy
- h5py
- Jupyter

Create the environment using:

```bash
git clone https://github.com/Griffin-Group/Angularly_Resolved_Modeling_of_Phonon_Interfacial_Transport.git
cd Angularly_Resolved_Modeling_of_Phonon_Interfacial_Transport
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

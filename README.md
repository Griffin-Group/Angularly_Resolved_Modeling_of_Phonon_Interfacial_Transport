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

1. **Acoustic Mismatch Model (AMM)** — Phonons are treated as plane waves that refract specularly at interfaces according to Snell’s law, \(\sin\theta_1/c_1 = \sin\theta_2/c_2\). The interface energy transmission is determined by the acoustic-impedance contrast:

$$\[
\alpha = 1-R,
\qquad
R =
\left(
\frac{Z_1\cos\theta_1-Z_2\cos\theta_2}
     {Z_1\cos\theta_1+Z_2\cos\theta_2}
\right)^2,
\qquad
Z=\rho c.
\]$$

This formulation captures the angular and polarization dependence shown in Figs. 5, 6, and S1.

2. **Tamura mass-difference scattering** — Random mass substitution in the SiGe alloy scatters phonons at a rate

$$\[
\tau_{\mathrm{md}}^{-1}(\omega)\propto \omega^2D(\omega),
\]$$

resulting in a frequency-dependent mean free path,

$$\[
\Lambda(\omega)=v_g(\omega)\tau_{\mathrm{md}}(\omega).
\]$$

The scattering layer is then treated as a phonon radiative-transport problem—the phonon analogue of the radiative transfer equation—and solved using the discrete-ordinates package **PythonicDISORT** (D. J. X. Ho; `pip install PythonicDISORT`) to obtain the transmittance spectrum shown in Figs. 2, 3, 4, S6, S7, and S8.

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

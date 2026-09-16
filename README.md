# Code and Data for “Theoretical Design of Phonon Filters for Quantum Sensing and Qubits”

This repository contains the code, input data, calculated results, and
figure-generation materials supporting the manuscript:

> **Theoretical Design of Phonon Filters for Quantum Sensing and Qubits**

**Manuscript status:** In preparation

The computational workflow implements an angularly resolved Acoustic Mismatch
Model (AMM) for phonon transmission and interfacial thermal transport. It is
used to calculate branch- and angle-resolved transmission probabilities,
effective transmission coefficients, and Kapitza conductance for material
interfaces relevant to superconducting quantum devices.

This repository is a pre-submission research release. The contents may be
updated before or during peer review.


## Physics

The AMM assumes specular (coherent) phonon scattering at a planar interface. For each acoustic branch i (TA1, TA2, LA), the transmission coefficient is given by:

$$\alpha_i = 1 - \left(\frac{Z_1\cos\theta_1 - Z_2\cos\theta_2}{Z_1\cos\theta_1 + Z_2\cos\theta_2}\right)^2$$

where Z = ρv is the acoustic impedance and θ₂ is the refraction angle from Snell's law (v₂/v₁ = sinθ₂/sinθ₁). Total internal reflection is handled automatically when sinθ₂ > 1. The interfacial heat flux is then:

$$q = \frac{1}{2}\sum_i v_i \int_0^{\omega_c}\int_0^1 \hbar\omega\, D_i(\omega)\,\alpha_i(\omega,\mu)\,\frac{\partial f}{\partial T}\,\mu\,d\mu\,d\omega$$

giving q in units of W m⁻² K⁻¹ (Kapitza conductance).

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

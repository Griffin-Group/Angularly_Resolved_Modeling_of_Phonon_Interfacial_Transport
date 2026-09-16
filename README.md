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


# Physics

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

# Repository Overview

1. **`AMM_jupyter-notebooks`**: Contains the Acoustic Mismatch Model (AMM) code, including an example calculation of the interfacial heat flux \(q\) between Si and Al.

2. **`Group_velocites_LA_TA`**: Contains the raw group-velocity data reported in this work, along with linear fits and their slopes for each material and propagation direction studied.

3. **`Inputs_outputs_data`**: Contains the DFT input and output files for all investigated materials, including the data used to plot the phonon dispersions.

4. **`Tamura_mass_difference_scattering`**: Demonstrates how DFT data are used as input to the Tamura mass-difference scattering model for the SiGe layer and how the frequency-dependent transmittance is calculated.

5. **`Transmission_coefficients_alpha`**: Demonstrates how the AMM is used to calculate angle- and polarization-dependent transmission coefficients in both the forward and reverse directions.

6. **`Figures`**: Contains all figures presented in the main text and Supplementary Information (SI) of the associated manuscript.

# Installation

## Requirements

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

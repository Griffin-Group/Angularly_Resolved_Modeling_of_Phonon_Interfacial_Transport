Here is a README summary:

---

# Acoustic Mismatch Model (AMM) for Phonon Transmission

## Overview

This repository implements the Acoustic Mismatch Model (AMM) to calculate phonon transmission coefficients and interfacial heat transfer across material interfaces relevant to superconducting qubit devices (e.g., Nb/SiGe, NbN/SiGe, Al/SiGe). The code computes angular-resolved transmission probabilities α(θ), effective transmission T_eff, and Kapitza conductance q (W/m²/K) from first-principles phonon data.

---

## Physics

The AMM assumes specular (coherent) phonon scattering at a planar interface. For each acoustic branch i (TA1, TA2, LA), the transmission coefficient is given by:

$$\alpha_i = 1 - \left(\frac{Z_1\cos\theta_1 - Z_2\cos\theta_2}{Z_1\cos\theta_1 + Z_2\cos\theta_2}\right)^2$$

where Z = ρv is the acoustic impedance and θ₂ is the refraction angle from Snell's law (v₂/v₁ = sinθ₂/sinθ₁). Total internal reflection is handled automatically when sinθ₂ > 1. The interfacial heat flux is then:

$$q = \frac{1}{2}\sum_i v_i \int_0^{\omega_c}\int_0^1 \hbar\omega\, D_i(\omega)\,\alpha_i(\omega,\mu)\,\frac{\partial f}{\partial T}\,\mu\,d\mu\,d\omega$$

giving q in units of W m⁻² K⁻¹ (Kapitza conductance).

---

## Code Structure

### `ThermalProperties` class
Loads phonon data for a single material from phonopy output (VASP or QE) and computes:
- **Density** ρ (kg/m³) — from primitive cell volume and atomic masses, with automatic Bohr→Å unit conversion for QE calculations (detected from `physical_unit: length: "au"` in `phonopy.yaml`)
- **Sound speed** v (m/s) — group velocity magnitude near Γ, averaged over 3 BZ directions for each acoustic branch
- **Mode-projected DOS** D(ω) (s/m³) — phonon density of states per branch, with Debye T³ correction at low frequencies
- **∂f/∂T** — derivative of Bose-Einstein distribution at the specified temperature (default 50 mK)

### `AMM` class
Combines two `ThermalProperties` objects (material_1 = source, material_2 = receiver) and computes:
- `get_alpha(θ)` — branch-resolved transmission coefficient vs incidence angle
- `get_effective_transmission()` — angle-averaged T_eff = 2∫α(μ)μ dμ
- `heat_transfer_coeff()` — Kapitza conductance q (W/m²/K)

### Standalone script (`simple_amm_alpha.py`)
Lightweight version requiring no phonopy — takes density and TA/LA group velocities directly as inputs. Outputs α(cosΘ) plot, T_eff, and saves data to a text file. Units can be m/s or km/s as long as they are consistent for both materials.

---

## Input Requirements

**Full AMM (phonopy-based):**
- `phonopy.yaml` — crystal structure, symmetry, supercell/primitive matrix, unit metadata
- `force_constants.hdf5` or `FORCE_SETS` — interatomic force constants from DFT (VASP or QE)
- `POSCAR` — unit cell structure

**Standalone script:**
- ρ₁, ρ₂ (kg/m³)
- v_TA, v_LA for each material (m/s or km/s, consistent)

---

## Unit Notes

| Quantity | Unit |
|---|---|
| Sound speed | m/s |
| Density | kg/m³ |
| Acoustic impedance Z = ρv | kg m⁻² s⁻¹ |
| Mode DOS D(ω) | s/m³ |
| Transmission α | dimensionless |
| Kapitza conductance q | W m⁻² K⁻¹ |
| T_eff | dimensionless |

**Important:** QE phonopy calculations use Bohr (au) for lengths. The code auto-detects this from `physical_unit: length: "au"` in `phonopy.yaml` and applies the correct Bohr³→Å³ conversion for volume-dependent quantities (density, DOS). Group velocities and frequencies are handled correctly by phonopy's internal unit conversion factor (`frequency_unit_conversion_factor` in the yaml).

---

## Supported Calculators

- **VASP** — lengths in Å, force constants in eV/Å²
- **QE (Quantum ESPRESSO)** — lengths in Bohr (au), force constants in Ry/Bohr², handled automatically

---

## Key Physical Notes

- Material 1 is always the **phonon source** — its DOS and sound speed enter the heat flux integral (Eq. 3). Material 2 contributes only through its acoustic impedance Z₂ in α.
- The TA and LA critical angles are in **opposite directions** when the velocity ordering flips between polarizations (e.g., v_LA(SiGe) > v_LA(NbN) but v_TA(SiGe) < v_TA(NbN)). This is physical, not a bug.
- Transmission is **not symmetric**: q(1→2) ≠ q(2→1) because the source DOS differs. Run both directions explicitly.
- Temperature is set to **50 mK** by default, relevant for superconducting qubit operating conditions.

#!/usr/bin/env python
"""
Convert QE q2r.x force constants to a phonopy-readable
phonopy_params_q2r.yaml file.

Fixes over the documented phonopy example script:
  1. FREQUENCY FACTOR: builds Phonopy with set_factor_by_calculator=True
     (or explicit ev-ang conversion) so the correct THz factor is used.
  2. CELL UNITS: read_pwscf() returns Bohr, not Angstrom. This script
     reports both units and optionally checks against expected values.
  3. MONOCLINIC / GENERAL CELLS: --expected-a/b/c checks each lattice
     vector independently, so non-cubic cells are validated correctly.

Usage:
    python make_fc_q2r.py scf.in nbn.fc [options]

Examples:
    # cubic P-43m NbN
    python make_fc_q2r.py scf.in nbn.fc --expected-a 4.422387 --units ev-ang

    # monoclinic C2/m NbN
    python make_fc_q2r.py scf.in nbn.fc \\
        --expected-a 11.199961 --expected-b 5.601403 --expected-c 11.201673 \\
        --units ev-ang
"""

import sys
import argparse
import numpy as np

from phonopy import Phonopy
from phonopy.interface.qe import read_pwscf, PH_Q2R
from phonopy.physical_units import get_calculator_physical_units, get_physical_units
from phonopy.structure.atoms import PhonopyAtoms


def check_lattice_param(name, parsed_ang, expected, tol, bohr_to_A):
    """Check one lattice constant (in Angstrom). Returns False and prints
    a diagnostic if the relative error exceeds tol."""
    if expected is None:
        return True
    rel_err = abs(parsed_ang - expected) / expected
    if rel_err > tol:
        ratio = parsed_ang / expected
        print(f"\nERROR: parsed |{name}| = {parsed_ang:.6f} Ang differs from "
              f"expected {expected:.6f} Ang by {rel_err*100:.2f}% "
              f"(tolerance {tol*100:.1f}%).")
        print(f"  ratio (parsed/expected) = {ratio:.6f}")
        if abs(ratio - 1.0 / bohr_to_A) < 1e-3:
            print("  This ratio matches 1/Bohr (~1.8897): likely a unit "
                  "keyword mismatch in scf.in (CELL_PARAMETERS or celldm).")
        elif abs(ratio - bohr_to_A) < 1e-3:
            print("  This ratio matches Bohr (~0.5292): same kind of "
                  "unit mismatch, opposite direction.")
        return False
    print(f"  |{name}| check OK: {parsed_ang:.6f} Ang "
          f"(relative error {rel_err*100:.4f}%)")
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("primcell_filename",
                        help="QE scf.in used for the PH/q2r run")
    parser.add_argument("q2r_filename",
                        help="q2r.x output file (e.g. nbn.fc)")
    parser.add_argument("--expected-a", type=float, default=None,
                        help="Expected |a| in Angstrom (sanity check)")
    parser.add_argument("--expected-b", type=float, default=None,
                        help="Expected |b| in Angstrom (sanity check)")
    parser.add_argument("--expected-c", type=float, default=None,
                        help="Expected |c| in Angstrom (sanity check)")
    parser.add_argument("--tol", type=float, default=0.01,
                        help="Relative tolerance for lattice checks (default 1%%)")
    parser.add_argument("--output", default="phonopy_params_q2r.yaml",
                        help="Output yaml filename")
    parser.add_argument("--units", choices=["qe", "ev-ang"], default="qe",
                        help=("'qe': keep native QE units (Bohr, Ry/au^2), "
                              "factor=108.970772. "
                              "'ev-ang': convert to Angstrom + eV/Ang^2, "
                              "factor=15.633302."))
    args = parser.parse_args()

    # ------------------------------------------------------------------ #
    # Step 1: read structure
    # read_pwscf() always returns in QE-native Bohr (au), by design.
    # ------------------------------------------------------------------ #
    cell, _ = read_pwscf(args.primcell_filename)
    pu = get_physical_units()
    bohr_to_A = pu.Bohr  # 0.52917721...

    lat_bohr = cell.cell
    lat_ang  = lat_bohr * bohr_to_A
    a_bohr, b_bohr, c_bohr = [np.linalg.norm(v) for v in lat_bohr]
    a_ang,  b_ang,  c_ang  = [np.linalg.norm(v) for v in lat_ang]
    volume_ang = abs(np.linalg.det(lat_ang))

    # Extract monoclinic angle beta (angle between a and c vectors)
    cos_beta = np.dot(lat_ang[0], lat_ang[2]) / (a_ang * c_ang)
    beta_deg = np.degrees(np.arccos(np.clip(cos_beta, -1, 1)))

    print("=" * 60)
    print(f"Parsed cell from: {args.primcell_filename}")
    print(f"  (read_pwscf native unit is Bohr/au)")
    print(f"  |a| = {a_bohr:.6f} au  =  {a_ang:.6f} Ang")
    print(f"  |b| = {b_bohr:.6f} au  =  {b_ang:.6f} Ang")
    print(f"  |c| = {c_bohr:.6f} au  =  {c_ang:.6f} Ang")
    print(f"  beta = {beta_deg:.4f} deg")
    print(f"  volume = {volume_ang:.6f} Ang^3")
    print(f"  natoms (primitive) = {len(cell.symbols)}")
    print("=" * 60)

    # ------------------------------------------------------------------ #
    # Step 2: sanity checks (any combination of a, b, c)
    # ------------------------------------------------------------------ #
    checks = [
        check_lattice_param("a", a_ang, args.expected_a, args.tol, bohr_to_A),
        check_lattice_param("b", b_ang, args.expected_b, args.tol, bohr_to_A),
        check_lattice_param("c", c_ang, args.expected_c, args.tol, bohr_to_A),
    ]
    if not all(checks):
        print("\nAborting — fix scf.in and rerun.")
        sys.exit(1)
    if any(x is not None for x in [args.expected_a, args.expected_b, args.expected_c]):
        print()

    # ------------------------------------------------------------------ #
    # Step 3: parse q2r force constants
    # ------------------------------------------------------------------ #
    q2r = PH_Q2R(args.q2r_filename)
    q2r.run(cell)

    # ------------------------------------------------------------------ #
    # Step 4: build Phonopy object with correct units
    # ------------------------------------------------------------------ #
    if args.units == "ev-ang":
        fc_factor = pu.Rydberg / pu.Bohr**2  # Ry/au^2 -> eV/Ang^2, ~48.5868

        prim = q2r.primitive
        cell_ang = PhonopyAtoms(
            symbols=prim.symbols,
            cell=prim.cell * bohr_to_A,
            scaled_positions=prim.scaled_positions,
            masses=prim.masses,
        )
        fc_ev_ang = q2r.fc * fc_factor

        print(f"Converting cell:            Bohr -> Angstrom   (x{bohr_to_A:.6f})")
        print(f"Converting force constants: Ry/au^2 -> eV/Ang^2 (x{fc_factor:.6f})")

        ph = Phonopy(
            cell_ang,
            supercell_matrix=q2r.dimension,
            calculator=None,  # eV/Ang^2 convention -> factor 15.633302
        )
        ph.force_constants = fc_ev_ang

    else:  # "qe"
        ph = Phonopy(
            q2r.primitive,
            supercell_matrix=q2r.dimension,
            calculator="qe",
            set_factor_by_calculator=True,  # forces 108.970772, not 15.633302
        )
        ph.force_constants = q2r.fc

    # ------------------------------------------------------------------ #
    # Step 5: verify factor and save
    # ------------------------------------------------------------------ #
    qe_factor = get_calculator_physical_units("qe").factor
    used_factor = ph.unit_conversion_factor
    print(f"Frequency conversion factor: {used_factor:.6f}", end="")
    if args.units == "qe":
        ok = abs(used_factor - qe_factor) < 1e-3
        print(f"  ({'OK' if ok else 'WARNING: expected ' + str(qe_factor)})")
    else:
        print()

    ph.save(args.output)
    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()

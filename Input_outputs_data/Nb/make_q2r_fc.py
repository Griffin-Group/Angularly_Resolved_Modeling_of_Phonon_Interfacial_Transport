#!/usr/bin/env python
"""
Convert QE q2r.x force constants to a phonopy-readable FORCE_CONSTANTS /
phonopy_params_q2r.yaml file, with two fixes applied on top of the
documented phonopy example script:

  1. FREQUENCY FACTOR: PH_Q2R.save() builds Phonopy(..., calculator="qe")
     without set_factor_by_calculator=True, so in phonopy <= ~2.4x it
     silently keeps the VASP-style default factor (15.633302) instead of
     the correct QE Ry/au^2 -> THz factor (108.970772). We build the
     Phonopy object ourselves and pass set_factor_by_calculator=True
     explicitly so this can't silently default to the wrong value.

  2. CELL UNITS: phonopy's read_pwscf() returns the unit cell in the
     QE-native length unit (Bohr/au), not Angstrom. If scf.in's
     CELL_PARAMETERS/celldm unit keyword doesn't match what the numbers
     actually are, you get a structure scaled by the Bohr<->Angstrom
     ratio (~1.8897x) without any error or warning. This script prints
     the parsed lattice constant and volume up front, and (optionally)
     checks it against an expected value you supply, so a units mistake
     in scf.in is caught here instead of silently propagating into every
     downstream phonon/group-velocity/AMM calculation.

Usage:
    python make_fc_q2r.py scf.in nbn.fc [--expected-a A_IN_ANGSTROM] [--tol 0.01]

Example:
    python make_fc_q2r.py scf.in nbn.fc --expected-a 4.422386913
"""

import sys
import argparse
import numpy as np

from phonopy import Phonopy
from phonopy.interface.qe import read_pwscf, PH_Q2R
from phonopy.physical_units import get_calculator_physical_units, get_physical_units
from phonopy.structure.atoms import PhonopyAtoms


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("primcell_filename", help="QE scf.in used for the PH/q2r run")
    parser.add_argument("q2r_filename", help="q2r.x output file (e.g. nbn.fc)")
    parser.add_argument(
        "--expected-a",
        type=float,
        default=None,
        help="Expected lattice constant 'a' in Angstrom, for a sanity check. "
        "If the parsed cell differs by more than --tol (relative), the "
        "script aborts before writing anything.",
    )
    parser.add_argument(
        "--tol",
        type=float,
        default=0.01,
        help="Relative tolerance for the --expected-a check (default 1%%).",
    )
    parser.add_argument(
        "--output",
        default="phonopy_params_q2r.yaml",
        help="Output yaml filename (default: phonopy_params_q2r.yaml)",
    )
    parser.add_argument(
        "--units",
        choices=["qe", "ev-ang"],
        default="qe",
        help="'qe' (default): keep native QE units (Bohr, Ry/au^2), tagged "
        "calculator='qe', factor=108.970772. 'ev-ang': explicitly convert "
        "the cell to Angstrom and force constants to eV/Angstrom^2, and "
        "save with the standard default factor=15.633302.",
    )
    args = parser.parse_args()

    # --- Step 1: read structure ---
    # NOTE: read_pwscf() returns the cell in QE's native length unit, which
    # is Bohr (au), by design -- not Angstrom. This is normal/expected, not
    # a units mistake in scf.in. We report both the raw (Bohr) and converted
    # (Angstrom) values below so the sanity check compares like with like.
    cell, _ = read_pwscf(args.primcell_filename)
    pu = get_physical_units()
    bohr_to_A = pu.Bohr  # 0.5291772109

    lat_bohr = cell.cell
    lat_ang = lat_bohr * bohr_to_A
    a_bohr, b_bohr, c_bohr = (np.linalg.norm(v) for v in lat_bohr)
    a_ang, b_ang, c_ang = (np.linalg.norm(v) for v in lat_ang)
    volume_ang = abs(np.linalg.det(lat_ang))

    print("=" * 60)
    print(f"Parsed cell from: {args.primcell_filename}")
    print(f"  (read_pwscf native unit is Bohr/au)")
    print(f"  |a| = {a_bohr:.6f} au  =  {a_ang:.6f} Ang")
    print(f"  |b| = {b_bohr:.6f} au  =  {b_ang:.6f} Ang")
    print(f"  |c| = {c_bohr:.6f} au  =  {c_ang:.6f} Ang")
    print(f"  volume = {volume_ang:.6f} Ang^3")
    print(f"  natoms (primitive) = {len(cell.symbols)}")
    print("=" * 60)

    if args.expected_a is not None:
        # Compare against the Angstrom value, since --expected-a is given
        # in Angstrom.
        rel_err = abs(a_ang - args.expected_a) / args.expected_a
        if rel_err > args.tol:
            ratio = a_ang / args.expected_a
            msg = [
                f"\nERROR: parsed |a| = {a_ang:.6f} Ang differs from",
                f"expected {args.expected_a:.6f} Ang by {rel_err*100:.2f}% "
                f"(tolerance {args.tol*100:.1f}%).",
                f"ratio (parsed/expected) = {ratio:.6f}",
            ]
            if abs(ratio - 1.0 / bohr_to_A) < 1e-3:
                msg.append(
                    "This ratio matches 1/Bohr (~1.8897): your scf.in likely "
                    "has CELL_PARAMETERS or celldm(1) in the wrong unit "
                    "keyword (angstrom numbers tagged as something else, or "
                    "vice versa). Fix scf.in and rerun."
                )
            elif abs(ratio - bohr_to_A) < 1e-3:
                msg.append(
                    "This ratio matches Bohr (~0.5292): same kind of "
                    "unit mismatch, opposite direction. Fix scf.in and rerun."
                )
            print("\n".join(msg))
            sys.exit(1)
        else:
            print(f"Cell check OK (relative error {rel_err*100:.3f}%).\n")

    # --- Step 2: parse q2r force constants ---
    q2r = PH_Q2R(args.q2r_filename)
    q2r.run(cell)

    if args.units == "ev-ang":
        # Explicit conversion: cell read_pwscf() returns is in Bohr (QE
        # native length unit); force constants from q2r are in Ry/au^2.
        # Convert both explicitly so the saved file is genuinely in
        # eV/Angstrom^2 / Angstrom / AMU, not just relabeled.
        fc_factor = pu.Rydberg / pu.Bohr**2       # Ry/au^2 -> eV/Ang^2, ~48.5868

        prim = q2r.primitive
        cell_ang = PhonopyAtoms(
            symbols=prim.symbols,
            cell=prim.cell * bohr_to_A,
            scaled_positions=prim.scaled_positions,
            masses=prim.masses,
        )
        fc_ev_ang = q2r.fc * fc_factor

        print(f"Converting cell: Bohr -> Angstrom (x{bohr_to_A:.6f})")
        print(f"Converting force constants: Ry/au^2 -> eV/Ang^2 (x{fc_factor:.6f})")

        ph = Phonopy(
            cell_ang,
            supercell_matrix=q2r.dimension,
            calculator=None,  # default -> eV/Angstrom^2, Angstrom, factor 15.633302
        )
        ph.force_constants = fc_ev_ang
    else:
        # --- Step 3: build Phonopy object with the correct QE frequency factor ---
        ph = Phonopy(
            q2r.primitive,
            supercell_matrix=q2r.dimension,
            calculator="qe",
            set_factor_by_calculator=True,  # <-- the actual fix for the 6.97x bug
        )
        ph.force_constants = q2r.fc

    qe_units = get_calculator_physical_units("qe")
    print(f"Frequency conversion factor used: {ph.unit_conversion_factor:.6f}")
    if args.units == "qe":
        print(f"(expected QE value: {qe_units.factor:.6f})")

    ph.save(args.output)
    print(f"\nSaved: {args.output}")


if __name__ == "__main__":
    main()

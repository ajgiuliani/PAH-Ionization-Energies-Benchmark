#!/usr/bin/env python3
"""
Physical consistency check of the raw energy files of the PAH benchmark.

    python3 check_results.py

The files checked are those of DEPOSITED, below: the three files of data/. The
script takes no argument; to check a fresh campaign, add its path to that list.

Each is a tab-separated file of 11 columns:

    Molecule  Method  E_neutral_Rn  E_cation_Rc  E_cation_Rn
              E_dication_Rd  E_dication_Rc  E_dication_Rn
              EPZ_neutral  EPZ_cation  EPZ_dication

WHAT IS CHECKED
---------------
The three relations tested are variational bounds: a relaxed geometry is, by
construction, the one that minimises the energy of its own species, hence

    E_cation(R_c)   <= E_cation(R_n)
    E_dication(R_d) <= E_dication(R_c)
    E_dication(R_d) <= E_dication(R_n)


"""
import os
import sys

BASE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(BASE, os.pardir, "data")

DEPOSITED = ["results_semi-empirical.txt", "results_g-xTB.txt",
             "results_composite.txt"]

# Tolerance in eV. The production scripts write the energies without rounding,
# so this threshold does not cover a loss of precision on writing, but the
# convergence residue of the optimiser, whose criterion is a force of
# 5e-4 eV/Angstrom.
TOL = 1e-4

# (label, indices to subtract) on the vector
#   0 E_neutral_Rn  1 E_cation_Rc  2 E_cation_Rn
#   3 E_dication_Rd 4 E_dication_Rc 5 E_dication_Rn
BOUNDS = [("E_cat(R_c) - E_cat(R_n)", 1, 2),
          ("E_di(R_d) - E_di(R_c)",   3, 4),
          ("E_di(R_d) - E_di(R_n)",   3, 5)]


def load(path):
    """Read a results file, tolerating header and tabulation variants."""
    out = {}
    with open(path) as f:
        for line in f.read().strip().split("\n")[1:]:
            p = line.rstrip("\t").split("\t")
            if len(p) >= 11:
                out[(p[0], p[1])] = [float(x) for x in p[2:11]]
    return out


def check(path):
    """Check one file. Returns (n_rows, n_errors)."""
    print(f"\n### {os.path.basename(path)}")
    try:
        rows = load(path)
    except OSError as exc:
        print(f"  !! illisible : {exc}")
        return 0, 1

    if not rows:
        print("  !! aucune ligne exploitable (11 colonnes attendues)")
        return 0, 1

    errors = 0
    for (mol, meth), v in sorted(rows.items()):
        for label, i, j in BOUNDS:
            delta = v[i] - v[j]
            if delta > TOL:
                errors += 1
                print(f"  !! {mol:14}{meth:11}{label} = {delta:+.4f} eV")

    verdict = "OK" if errors == 0 else f"{errors} violation(s)"
    print(f"  {len(rows):3d} lignes controlees - {verdict}")
    return len(rows), errors


def main():
    paths = [os.path.join(DATA, f) for f in DEPOSITED]

    total = errors = 0
    for path in paths:
        n, e = check(path)
        total += n
        errors += e

    print(f"\n{total} lignes, {errors} violation(s) de borne variationnelle")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())

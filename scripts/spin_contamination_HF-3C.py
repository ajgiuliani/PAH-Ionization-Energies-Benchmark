#!/usr/bin/env python3
"""
Spin contamination diagnostic of the PAH cations in HF-3C.

For each PAH of the benchmark, runs single points on the HF-3C geometries
already optimised:

  - doublet cation at the geometry of the cation  (R_c): UHF then ROHF
  - doublet cation at the geometry of the neutral (R_n): UHF then ROHF
  - coronene only: triplet dication at R_d, R_c, R_n: UHF then ROHF

and extracts from each ORCA output:
  - <S**2> and its ideal value S(S+1)
  - the total energy

Produces a TSV file with, per species:
  S2_UHF, S2_ideal, contamination = S2_UHF - S2_ideal,
  E_UHF, E_ROHF, dE = E_ROHF - E_UHF   (variational stabilisation of the UHF
  symmetry breaking)

No geometry optimisation, single points only.

Usage:
    python3 spin_contamination_HF-3C.py
"""

import os
import re
import time

from ase.io import read
from ase.calculators.orca import ORCA, OrcaProfile

# ---------------------------------------------------------------- configuration

# ORCA is rarely on the PATH, and /usr/bin/orca is the GNOME screen reader on
# most distributions rather than the quantum chemistry program: the full path is
# mandatory. Set ORCA_BIN in the environment to point at the executable.
ORCA_BIN = os.environ.get("ORCA_BIN", "orca")
NPROCS = 32
MAXCORE = 4000
METHOD = "HF-3C"

STRUCT_DIR = "structures/optimized"

NAMES = ["naphthalene", "fluorene", "anthracene", "pyrene",
         "tetracene", "coronene", "ovalene", "dicoronylene"]

# The coronene dication has a triplet ground state, the others are singlets.
TRIPLET_DICATION = ["coronene"]

profile = OrcaProfile(command=ORCA_BIN)
os.environ["ASE_ORCA_COMMAND"] = ORCA_BIN + " PREFIX.inp > PREFIX.out"

BLOCK_UHF = f"%pal nprocs {NPROCS} end \n%maxcore {MAXCORE}"
BLOCK_ROHF = f"%scf HFTyp ROHF end \n%pal nprocs {NPROCS} end \n%maxcore {MAXCORE}"

# ---------------------------------------------------------------- extraction

RE_S2 = re.compile(r"Expectation value of <S\*\*2>\s*:\s*([-\d.]+)")
RE_S2_IDEAL = re.compile(r"Ideal value S\*\(S\+1\) for S=[\d.]+\s*:\s*([-\d.]+)")


def parse_s2(outfile):
    """Return (<S**2>, ideal value) of the last occurrence, or (None, None).

    BEWARE: with ASE 3.23 the OrcaTemplate imposes _label = 'orca', so the
    output is always called orca.out whatever the label passed to the
    calculator. That is why each single point runs in its own directory and why
    we read <directory>/orca.out.
    """
    try:
        txt = open(outfile, errors="ignore").read()
    except OSError:
        return None, None
    s2 = RE_S2.findall(txt)
    ideal = RE_S2_IDEAL.findall(txt)
    return (float(s2[-1]) if s2 else None,
            float(ideal[-1]) if ideal else None)


def single_point(atoms, charge, mult, workdir, block):
    """ORCA single point in its own directory.

    Returns (energy_eV, <S**2>, S2_ideal), or (None, None, None) on failure.
    The dedicated directory serves two purposes: keeping the output of each
    species (ASE always writes orca.out and would otherwise overwrite it), and
    being able to read it back.
    """
    os.makedirs(workdir, exist_ok=True)
    atoms = atoms.copy()
    atoms.center()
    atoms.calc = ORCA(charge=charge, mult=mult,
                      directory=workdir,
                      profile=profile,
                      orcasimpleinput=METHOD,
                      orcablocks=block)
    try:
        energy = atoms.get_potential_energy()
    except Exception as exc:                      # SCF not converged, ORCA crash...
        print(f"    ECHEC {workdir}: {exc}")
        return None, None, None
    s2, ideal = parse_s2(os.path.join(workdir, "orca.out"))
    return energy, s2, ideal


def geometry(name, kind):
    """kind among 'neutral', 'cation', 'dication'."""
    path = f"{STRUCT_DIR}/{name}_{METHOD}_{kind}_bfgs.xyz"
    if not os.path.isfile(path):
        print(f"    structure absente : {path}")
        return None
    return read(path)


# ----------------------------------------------------------------- main loop

def main():
    fname = "spin_contamination_" + time.strftime("%Y%m%d-%H%M%S") + ".txt"
    with open(fname, "w") as f:
        f.write("\t".join(["Molecule", "Species", "Geometry", "Mult",
                           "S2_UHF", "S2_ideal", "Contamination",
                           "E_UHF", "E_ROHF", "dE_ROHF_minus_UHF"]) + "\n")

    print(f"Sortie : {fname}\n")

    for name in NAMES:
        print(f"--- {name}")

        # species to treat: (label, geometry, charge, multiplicity)
        jobs = [("cation", "cation", 1, 2),
                ("cation", "neutral", 1, 2)]
        if name in TRIPLET_DICATION:
            jobs += [("dication", "dication", 2, 3),
                     ("dication", "cation", 2, 3),
                     ("dication", "neutral", 2, 3)]

        for species, geom_kind, charge, mult in jobs:
            atoms = geometry(name, geom_kind)
            if atoms is None:
                continue

            tag = f"{name}_{species}_R{geom_kind[0]}"
            print(f"  {species} @ R_{geom_kind[0]} (mult={mult})")

            e_uhf, s2, s2_ideal = single_point(atoms, charge, mult,
                                               f"spin/{tag}_uhf", BLOCK_UHF)
            e_rohf, _, _ = single_point(atoms, charge, mult,
                                        f"spin/{tag}_rohf", BLOCK_ROHF)

            contamination = (s2 - s2_ideal) if (s2 is not None and s2_ideal is not None) else None
            delta = (e_rohf - e_uhf) if (e_uhf is not None and e_rohf is not None) else None

            def fmt(x):
                return f"{x:.6f}" if x is not None else "NA"

            row = [name, species, f"R_{geom_kind[0]}", str(mult),
                   fmt(s2), fmt(s2_ideal), fmt(contamination),
                   fmt(e_uhf), fmt(e_rohf), fmt(delta)]
            with open(fname, "a") as f:
                f.write("\t".join(row) + "\n")

            print(f"    <S**2> = {fmt(s2)} (ideal {fmt(s2_ideal)}), "
                  f"E_ROHF - E_UHF = {fmt(delta)} eV")

    print(f"\nTermine. Resultats dans {fname}")


if __name__ == "__main__":
    main()

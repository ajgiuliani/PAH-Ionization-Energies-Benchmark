#!/usr/bin/env python3
"""
g-xTB production campaign for the ionization energies of PAHs.

For each molecule, the six energies from which the adiabatic and vertical
ionization energies follow are computed:

    E_neutral_Rn    neutral  at its own geometry            (R_n)
    E_cation_Rc     cation   at its own geometry            (R_c)
    E_cation_Rn     cation   at the geometry of the neutral (R_n)
    E_dication_Rd   dication at its own geometry            (R_d)
    E_dication_Rc   dication at the geometry of the cation  (R_c)
    E_dication_Rn   dication at the geometry of the neutral (R_n)

The three species are relaxed one after the other, each starting from the
converged geometry of the previous one: the cation starts from the relaxed
neutral, the dication from the relaxed cation. This is the cheapest warm start,
the geometries of the three species remaining close.

NO ZERO-POINT ENERGY
--------------------
The three EPZ columns are written as 0.0

SPIN STATE
----------
The number of unpaired electrons is passed through the uhf keyword, which the
calculator writes into the .UHF file read by g-xTB: 0 for the neutral, 1 for the
cation, 0 for the singlet dications and 2 for that of coronene, whose ground
state is a triplet.

DEPENDENCIES
------------
The ASE calculator for g-xTB, gxtb.py, is provided next to this script and must
be installed as ase/calculators/gxtb.py for the import to work. It builds its
command line from the GXTB_COMMAND variable and writes the charge and the
multiplicity into the .CHRG and .UHF files.

Usage:
    python3 benchmark-gxtb.py

to be run from the directory holding structures/initial/, where the starting
geometries of the astrochemical database are found.
"""

import os
import time

from ase.io import read, write
from ase.optimize import BFGS
from ase.calculators.gxtb import Gxtb

# Location of the g-xTB installation. GXTBHOME points to the parameter files,
# GXTB_COMMAND to the binary; the calculator adds the -c coord -grad options
# itself.
#
# The tilde is expanded here because neither g-xTB nor the shell that calls it
# would do so. The two variables are only set if they are not set already: an
# installation elsewhere is declared in the environment, without touching the
# script.
GXTB_HOME = os.path.expanduser('~/g-xtb')
os.environ.setdefault('GXTBHOME', f'{GXTB_HOME}/parameters/')
os.environ.setdefault('GXTB_COMMAND', f'{GXTB_HOME}/binary/gxtb')


# ---------------------------------------------------------------- parameters

MOLECULES = ['naphthalene', 'fluorene', 'anthracene', 'pyrene',
             'tetracene', 'coronene', 'ovalene', 'dicoronylene']

METHOD = 'g-xTB'

# Molecules whose ground-state dication is a triplet. The assignment is taken
# from the database of Malloci et al.; for the seven other molecules the
# dication is a closed-shell singlet.
TRIPLET_DICATION = ['coronene']

MAX_FORCE = 0.0005                  # BFGS convergence criterion, eV/Angstrom

INITIAL_DIR = 'structures/initial'
OPTIMIZED_DIR = 'structures/optimized'

# Restart file written by g-xTB in the current directory. It is deleted after
# each molecule: kept from one molecule to the next, it would restart the SCF
# from a density that no longer matches the geometry read.
RESTART_FILE = 'gxtbrestart'

CHARGE = {'neutral': 0, 'cation': 1, 'dication': 2}


# --------------------------------------------------------------- input-output

def optimized_path(molecule, species, ext='traj'):
    return f'{OPTIMIZED_DIR}/{molecule}_{METHOD}_{species}_bfgs.{ext}'


def load_start(molecule, species, fallback=None):
    """Starting geometry of a species.

    Restart after an interruption: if the trajectory of a previous optimisation
    exists, we resume from its last image and BFGS stops at the first step.
    Otherwise we start from fallback, the relaxed geometry of the previous
    species, or from the geometry of the astrochemical database for the neutral.

    The fallback is copied: the optimisation of this species must not move the
    atoms of the previous one, whose energies at the non-relaxed geometries
    remain to be computed below.
    """
    path = optimized_path(molecule, species)
    if os.path.isfile(path):
        print('  geometrie optimisee trouvee :', path, flush=True)
        return read(path)
    if fallback is not None:
        print('  depart de la geometrie relaxee precedente', flush=True)
        return fallback.copy()
    path = f'{INITIAL_DIR}/{molecule}.xyz'
    print('  geometrie initiale           :', path, flush=True)
    return read(path)


def saveenergies(fname, molecule, energies, zpe):
    """Append one line to the results file, in the order of the header.

    The energies are in eV, as returned by ASE. We write str(x) without
    rounding: the ionization differences bear on the last decimals of numbers
    of the order of 10^4 eV, and any rounding here would be destructive.
    """
    values = [energies[k] for k in ('E_neutral_Rn', 'E_cation_Rc',
                                    'E_cation_Rn', 'E_dication_Rd',
                                    'E_dication_Rc', 'E_dication_Rn')]
    values += [zpe[k] for k in ('neutral', 'cation', 'dication')]
    with open(fname, 'a') as f:
        # no trailing tabulation: it would add an empty column
        f.write('\t'.join([molecule, METHOD] + [str(v) for v in values]) + '\n')


def write_header(fname):
    columns = ['Molecule', 'Method',
               'E_neutral_Rn', 'E_cation_Rc', 'E_cation_Rn',
               'E_dication_Rd', 'E_dication_Rc', 'E_dication_Rn',
               'EPZ_neutral', 'EPZ_cation', 'EPZ_dication']
    with open(fname, 'w') as f:
        f.write('\t'.join(columns) + '\n')


# ------------------------------------------------------------------ the calculation

def relax(atoms, molecule, species, uhf):
    """Relax a species, save its geometry and return its energy."""
    print(f'\n  optimisation du {species}', flush=True)
    atoms.center()
    atoms.calc = Gxtb(charge=CHARGE[species], uhf=uhf)
    dyn = BFGS(atoms, trajectory=optimized_path(molecule, species))
    dyn.run(fmax=MAX_FORCE)
    write(optimized_path(molecule, species, ext='xyz'), atoms)
    return atoms.get_potential_energy()


def single_point(atoms, species, uhf):
    """Energy of a given charge at the geometry of the object passed.

    A fresh calculator is built at each call, which guarantees that no result is
    recycled from one point to the next.
    """
    atoms.calc = Gxtb(charge=CHARGE[species], uhf=uhf)
    return atoms.get_potential_energy()


def run(molecule, fname):
    uhf_dication = 2 if molecule in TRIPLET_DICATION else 0
    if uhf_dication:
        print('         dication traite en triplet', flush=True)

    energies = {}

    # The three relaxed species, each starting from the previous one.
    neutral = load_start(molecule, 'neutral')
    energies['E_neutral_Rn'] = relax(neutral, molecule, 'neutral', 0)

    cation = load_start(molecule, 'cation', fallback=neutral)
    energies['E_cation_Rc'] = relax(cation, molecule, 'cation', 1)

    dication = load_start(molecule, 'dication', fallback=cation)
    energies['E_dication_Rd'] = relax(dication, molecule, 'dication',
                                      uhf_dication)

    # The three energies at the non-relaxed geometries: we change the charge of
    # the calculator without touching the positions.
    energies['E_cation_Rn'] = single_point(neutral, 'cation', 1)
    energies['E_dication_Rc'] = single_point(cation, 'dication', uhf_dication)
    energies['E_dication_Rn'] = single_point(neutral, 'dication', uhf_dication)

    for label in ('E_neutral_Rn', 'E_cation_Rc', 'E_cation_Rn',
                  'E_dication_Rd', 'E_dication_Rc', 'E_dication_Rn'):
        print(f'  {label:16} {energies[label]:18.9f} eV', flush=True)

    # no frequencies: see the header of this file
    zpe = {'neutral': 0.0, 'cation': 0.0, 'dication': 0.0}

    saveenergies(fname, molecule, energies, zpe)

    if os.path.isfile(RESTART_FILE):
        os.remove(RESTART_FILE)


# ----------------------------------------------------------------- main loop

def main():
    os.makedirs(OPTIMIZED_DIR, exist_ok=True)

    fname = 'results_' + time.strftime('%Y%m%d-%H%M%S') + '.txt'
    write_header(fname)
    print('resultats dans', fname, flush=True)

    print('\n' + '-' * 17 + METHOD + '-' * 17, flush=True)
    for molecule in MOLECULES:
        print('\n         ' + molecule, flush=True)
        run(molecule, fname)
        print('\n', flush=True)

    print('\n\n -------- FIN DES CALCULS -------- \n')


if __name__ == '__main__':
    main()

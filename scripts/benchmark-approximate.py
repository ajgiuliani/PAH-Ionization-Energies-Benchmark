#!/usr/bin/env python3
"""
Production run of the approximate methods for the ionization energies of PAHs.

THE SIX ENERGIES
----------------
For each molecule and each method:

    E_neutral_Rn    neutral  at its own geometry            (R_n)
    E_cation_Rc     cation   at its own geometry            (R_c)
    E_cation_Rn     cation   at the geometry of the neutral (R_n)
    E_dication_Rd   dication at its own geometry            (R_d)
    E_dication_Rc   dication at the geometry of the cation  (R_c)
    E_dication_Rn   dication at the geometry of the neutral (R_n)

plus the three zero-point energies, obtained by finite differences of the
forces (ase.vibrations, 6N + 1 gradient evaluations per species).

The three species are relaxed; the three energies at the non-relaxed geometries
are single points obtained by changing the charge of the calculator without
touching the positions.

The adiabatic and vertical ionization energies follow from these six numbers,
with the zero-point corrections where they make sense; that calculation is done
downstream by analysis.jl.

THE METHODS
-----------
  DFTB2   DFTB+ self-consistent charge, Slater-Koster set mio-1-1
  DFTB3   the same with full third order, set 3ob-3-1, Hubbard derivatives
  IPEA1   IPEA1-xTB, the xTB Hamiltonian of DFTB+
  PM6     \
  PM7      >  MOPAC
  RM1     /

All the DFTB+ methods run with a Fermi filling at 300 K and 250 SCC iterations
at most.

SPIN STATE OF THE DICATIONS
---------------------------
The dications are treated closed-shell (singlet), except that of coronene,
whose ground state is a triplet. For that molecule alone, the three single
points of charge 2 and the optimisation of the dication go open-shell:
colinear spin polarisation with two unpaired electrons on the DFTB+ side,
TRIPLET UHF keyword on the MOPAC side.

The cations, for their part, have an odd number of electrons but are not spin
polarised: DFTB+ fills the frontier orbital by half and MOPAC applies its
half-electron method. This is the usual restricted treatment at this level of
theory, and it is the one used for the article; it should be kept in mind when
comparing with the composite methods, where the cation is explicitly UHF.

DEPENDENCY
----------
The mult keyword of MOPAC('TRIPLET UHF') does not exist in the distributed ASE:
it comes from a local addition to ase/calculators/mopac.py, which adds the ms
and mult parameters and forces the DISP keyword so that the total energy is
available. Without that fix, the coronene dication would be computed
closed-shell and the keyword would be silently ignored.

Usage:
    python3 benchmark-approximate.py

to be run from the directory holding structures/initial/, where the starting
geometries of the astrochemical database are found.
"""

import os
import time

from ase.io import read, write
from ase.optimize import BFGS
from ase.vibrations import Vibrations
from ase.calculators.dftb import Dftb
from ase.calculators.mopac import MOPAC


# ---------------------------------------------------------------- parameters

MOLECULES = ['naphthalene', 'fluorene', 'anthracene', 'pyrene',
             'tetracene', 'coronene', 'ovalene', 'dicoronylene']

METHODS = ['DFTB2', 'DFTB3', 'IPEA1', 'PM6', 'PM7', 'RM1']

MOPAC_METHODS = ['PM6', 'PM7', 'RM1']

# Molecules whose ground-state dication is a triplet. The assignment is taken
# from the database of Malloci et al.; for the seven other molecules the
# dication is a closed-shell singlet.
TRIPLET_DICATION = ['coronene']

MAX_FORCE = 0.0005                  # BFGS convergence criterion, eV/Angstrom

# Directories of the DFTB+ Slater-Koster sets. The tilde is expanded here: ASE
# copies this path as it stands into dftb_in.hsd and DFTB+ would not interpret
# it. An installation elsewhere is declared through DFTB_SK_HOME, without
# touching the script.
SK_HOME = os.path.expanduser(os.environ.get('DFTB_SK_HOME', '~/dftb+'))
SK_PREFIX = {
    'DFTB2': f'{SK_HOME}/mio-1-1/',
    'DFTB3': f'{SK_HOME}/3ob-3-1/',
}

# Hubbard derivatives of the 3ob set, indispensable to the third order.
HUBBARD_DERIVS = {'C': '-0.1492', 'H': '-0.1857'}

# xTB Hamiltonian of DFTB+: the name of the parametrisation, not a set of .skf
# files.
XTB_METHOD = {'IPEA1': 'IPEA1-xTB'}

# Atomic spin constants, one value per atom (ShellResolvedSpin = No), proper to
# each parametrisation. They only come into play for the spin-polarised
# species, hence here the triplet dication alone.
SPIN_CONSTANTS = {
    'DFTB2': {'C': '-0.022',    'H': '-0.064'},
    'DFTB3': {'C': '-0.02265',  'H': '-0.07174'},
    'IPEA1': {'C': '-0.022725', 'H': '-0.071550'},
}

# Electronic temperature of 300 K expressed in Hartree (kB = 3.16681e-6 Ha/K).
# The Fermi filling it opens serves to unblock the SCC of the dications, whose
# frontier orbitals are nearly degenerate.
ELECTRONIC_TEMPERATURE = '0.950045E-03'
MAX_SCC_ITERATIONS = '250'

INITIAL_DIR = 'structures/initial'
OPTIMIZED_DIR = 'structures/optimized'
VIB_DIR = 'vib'

# Suffix of the initial geometry file of each species, in the nomenclature of
# the astrochemical database. The triplet dication of coronene also starts from
# the singlet geometry, the only starting geometry available for every
# molecule; it is no more than a starting point for BFGS.
INITIAL_SUFFIX = {'neutral': '', 'cation': '+', 'dication': '++s'}

CHARGE = {'neutral': 0, 'cation': 1, 'dication': 2}


# ------------------------------------------------------------------- factory

def calculator(method, charge, unpaired=0):
    """Return the ASE calculator of a method, for a given charge.

    charge   : 0, 1 or 2 electrons removed.
    unpaired : number of unpaired electrons. 0 asks for a closed-shell
               treatment or, for an odd number of electrons, for the default
               restricted treatment of the code; 2 asks for an explicitly
               spin-polarised triplet.

    One single function for the six methods: it is the only way to see how they
    really differ.
    """
    if method in MOPAC_METHODS:
        # MOPAC only writes the CHARGE keyword if the charge is non-zero, and
        # handles the odd number of electrons of the cation by itself.
        kwargs = {'method': method, 'charge': float(charge)}
        if unpaired == 2:
            # mult parameter added locally to ase/calculators/mopac.py.
            kwargs['mult'] = 'TRIPLET UHF'
        return MOPAC(**kwargs)

    # DFTB+ family. The keys are written into the dftb_in.hsd file in
    # alphabetical order by ASE, so the order of the dictionary is of no
    # importance; the trailing underscore marks the opening of a block.
    par = {
        'Hamiltonian_MaxAngularMomentum_': '',
        'Hamiltonian_MaxAngularMomentum_C': '"p"',
        'Hamiltonian_MaxAngularMomentum_H': '"s"',
        'Hamiltonian_Charge': f'{float(charge):.1f}',
        'Hamiltonian_MaxSCCIterations': MAX_SCC_ITERATIONS,
        'Hamiltonian_Filling_': 'Fermi',
        'Hamiltonian_Filling_Temperature': ELECTRONIC_TEMPERATURE,
        'ParserOptions_': '',
        'ParserOptions_ParserVersion': '10',
        'Parallel': '{ UseOmpThreads = Yes }',
    }

    if method in XTB_METHOD:
        # xTB defines its own basis: the MaxAngularMomentum block above makes
        # no sense here and DFTB+ ignores it, IgnoreUnprocessedNodes being Yes.
        # It is kept as it stands, the published campaign having passed it too.
        par['Hamiltonian_'] = 'xTB'
        par['Hamiltonian_Method'] = XTB_METHOD[method]
    else:
        par['Hamiltonian_'] = 'DFTB'
        par['Hamiltonian_SCC'] = 'Yes'
        par['Hamiltonian_SlaterKosterFiles_'] = 'Type2FileNames'
        par['Hamiltonian_SlaterKosterFiles_Prefix'] = SK_PREFIX[method]
        par['Hamiltonian_SlaterKosterFiles_Separator'] = '"-"'
        par['Hamiltonian_SlaterKosterFiles_Suffix'] = '".skf"'
        if method == 'DFTB3':
            par['Hamiltonian_ThirdOrderFull'] = 'Yes'
            par['Hamiltonian_HubbardDerivs_'] = ''
            for symbol, value in HUBBARD_DERIVS.items():
                par[f'Hamiltonian_HubbardDerivs_{symbol}'] = value

    if unpaired:
        par['Hamiltonian_SpinPolarisation_'] = 'Colinear'
        par['Hamiltonian_SpinPolarisation_UnpairedElectrons'] = str(unpaired)
        par['Hamiltonian_SpinConstants_'] = ''
        for symbol, value in SPIN_CONSTANTS[method].items():
            par[f'Hamiltonian_SpinConstants_{symbol}'] = value
        # Default value of DFTB+, written explicitly to remove the ambiguity:
        # one constant per atom and not per shell.
        par['Hamiltonian_SpinConstants_ShellResolvedSpin'] = 'No'

    return Dftb(**par)


# --------------------------------------------------------------- input-output

def optimized_path(molecule, method, species, ext='traj'):
    return f'{OPTIMIZED_DIR}/{molecule}_{method}_{species}_bfgs.{ext}'


def load_initial(molecule, method, species):
    """Starting geometry of a species.

    Restart after an interruption: if the trajectory of a previous optimisation
    exists, we resume from its last image and BFGS stops at the first step.
    Otherwise we start from the geometry of the astrochemical database.
    """
    path = optimized_path(molecule, method, species)
    if os.path.isfile(path):
        print('  geometrie optimisee trouvee :', path)
        return read(path)
    path = f'{INITIAL_DIR}/{molecule}{INITIAL_SUFFIX[species]}.xyz'
    print('  geometrie initiale           :', path)
    return read(path)


def load_optimized(molecule, method, species):
    """Converged geometry of a species, as it served for the energies."""
    return read(optimized_path(molecule, method, species))


def saveenergies(fname, molecule, method, energies, zpe):
    """Append one line to the results file, in the order of the header.

    The energies are in eV, as returned by ASE. We write str(x) without
    rounding: the ionization differences bear on the last decimals of numbers
    of the order of 10^2 to 10^4 eV, and any rounding here would be
    destructive.
    """
    values = [energies[k] for k in ('E_neutral_Rn', 'E_cation_Rc',
                                    'E_cation_Rn', 'E_dication_Rd',
                                    'E_dication_Rc', 'E_dication_Rn')]
    values += [zpe[k] for k in ('neutral', 'cation', 'dication')]
    with open(fname, 'a') as f:
        # no trailing tabulation: it would add an empty column
        f.write('\t'.join([molecule, method] + [str(v) for v in values]) + '\n')


def write_header(fname):
    columns = ['Molecule', 'Method',
               'E_neutral_Rn', 'E_cation_Rc', 'E_cation_Rn',
               'E_dication_Rd', 'E_dication_Rc', 'E_dication_Rn',
               'EPZ_neutral', 'EPZ_cation', 'EPZ_dication']
    with open(fname, 'w') as f:
        f.write('\t'.join(columns) + '\n')


# ------------------------------------------------------------ the three steps

def optimize(molecule, method, unpaired_dication):
    """Relax the three species and save the converged geometries."""
    for species in ('neutral', 'cation', 'dication'):
        print(f'\n  optimisation du {species}')
        atoms = load_initial(molecule, method, species)
        atoms.center()
        unpaired = unpaired_dication if species == 'dication' else 0
        atoms.calc = calculator(method, CHARGE[species], unpaired)
        dyn = BFGS(atoms, trajectory=optimized_path(molecule, method, species))
        dyn.run(fmax=MAX_FORCE)
        write(optimized_path(molecule, method, species, ext='xyz'), atoms)


def zero_point_energies(molecule, method, unpaired_dication):
    """Zero-point energies of the three relaxed species.

    ase.vibrations caches the 6N + 1 gradients in a directory named after the
    molecule, the method and the species: rerunning the script does not redo
    the calculations already done, but a cache left by a different geometry or
    parametrisation would be read back as it stands. When in doubt, delete the
    corresponding vib/ directory.
    """
    zpe = {}
    for species in ('neutral', 'cation', 'dication'):
        print(f'\n  frequences du {species}')
        atoms = load_optimized(molecule, method, species)
        unpaired = unpaired_dication if species == 'dication' else 0
        atoms.calc = calculator(method, CHARGE[species], unpaired)
        vib = Vibrations(atoms, name=f'{VIB_DIR}/vib_{molecule}_{method}_{species}')
        vib.run()
        zpe[species] = vib.get_zero_point_energy()
    return zpe


def total_energies(molecule, method, unpaired_dication):
    """The six total energies.

    Each species is read back from its converged geometry, then the wanted
    charges are applied to it in turn: the neutral carries the charges 0, 1 and
    2, the cation the charges 1 and 2, the dication the charge 2. A fresh
    calculator is built each time, which guarantees that no result is recycled
    from one point to the next.
    """
    # (species whose geometry is taken, charge applied, name of the column)
    points = [('neutral',  0, 'E_neutral_Rn'),
              ('cation',   1, 'E_cation_Rc'),
              ('neutral',  1, 'E_cation_Rn'),
              ('dication', 2, 'E_dication_Rd'),
              ('cation',   2, 'E_dication_Rc'),
              ('neutral',  2, 'E_dication_Rn')]

    energies = {}
    for species, charge, label in points:
        atoms = load_optimized(molecule, method, species)
        # the spin polarisation follows the charge, not the geometry: the three
        # points of charge 2 of coronene are triplets
        unpaired = unpaired_dication if charge == 2 else 0
        atoms.calc = calculator(method, charge, unpaired)
        energies[label] = atoms.get_potential_energy()
        print(f'  {label:16} {energies[label]:18.9f} eV')
    return energies


# ----------------------------------------------------------------- main loop

def main():
    os.makedirs(OPTIMIZED_DIR, exist_ok=True)
    os.makedirs(VIB_DIR, exist_ok=True)

    fname = 'results_' + time.strftime('%Y%m%d-%H%M%S') + '.txt'
    write_header(fname)
    print('resultats dans', fname)

    for method in METHODS:
        print('\n' + '-' * 17 + method + '-' * 17)

        for molecule in MOLECULES:
            print('\n         ' + molecule)
            unpaired_dication = 2 if molecule in TRIPLET_DICATION else 0
            if unpaired_dication:
                print('         dication traite en triplet')

            optimize(molecule, method, unpaired_dication)
            zpe = zero_point_energies(molecule, method, unpaired_dication)
            print('\n  EPZ neutre {neutral:.6f}  cation {cation:.6f}'
                  '  dication {dication:.6f}'.format(**zpe))

            energies = total_energies(molecule, method, unpaired_dication)
            saveenergies(fname, molecule, method, energies, zpe)

    print('\n\n -------- FIN DES CALCULS -------- \n')


if __name__ == '__main__':
    main()

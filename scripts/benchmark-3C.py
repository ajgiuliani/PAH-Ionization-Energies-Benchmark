#!/usr/bin/env python3
"""
Production run of the composite methods for the ionization energies of PAHs.

For each molecule and each method, it computes the six energies from which the
adiabatic and vertical ionization energies follow:

    E_neutral_Rn    neutral  at its own geometry            (R_n)
    E_cation_Rc     cation   at its own geometry            (R_c)
    E_cation_Rn     cation   at the geometry of the neutral (R_n)
    E_dication_Rd   dication at its own geometry            (R_d)
    E_dication_Rc   dication at the geometry of the cation  (R_c)
    E_dication_Rn   dication at the geometry of the neutral (R_n)

plus the three zero-point energies, obtained by finite differences of the forces
(ase.vibrations, that is 6N + 1 gradient evaluations per species).

The three species are relaxed; the three energies at the non-relaxed geometries
are single points obtained by changing the charge of the calculator without
touching the positions, which is why the neutral and cation objects are reused
further down in the loop.

Usage:
    python3 benchmark-3C.py
"""

import os
import time

from ase.io import read, write
from ase.optimize import BFGS
from ase.vibrations import Vibrations
from ase.calculators.orca import ORCA, OrcaProfile

# ORCA is rarely on the PATH, and /usr/bin/orca is the GNOME screen reader on
# most distributions rather than the quantum chemistry program: the full path is
# mandatory. Set ORCA_BIN in the environment to point at the executable.
ORCA_BIN = os.environ.get("ORCA_BIN", "orca")
profile = OrcaProfile(command=ORCA_BIN)
os.environ["ASE_ORCA_COMMAND"] = f"{ORCA_BIN} PREFIX.inp > PREFIX.out"

#---------Functions

def saveenergies(n, method, fname, E_neutral_Rn, E_cation_Rc, E_cation_Rn,
                 E_dication_Rd, E_dication_Rc, E_dication_Rn,
                 EPZ_neutral, EPZ_cation, EPZ_dication):
    """Append one line to the results file, in the order of the header.

    The energies are in eV, as returned by ASE. We write str(x) without
    rounding: the ionization differences bear on the sixth decimal of numbers
    of the order of 10^4 eV, and any rounding here would be destructive.
    """
    f = open(fname, "a")
    f.write(n+'\t')
    f.write(method+'\t')
    f.write(str( E_neutral_Rn  )+'\t')
    f.write(str( E_cation_Rc   )+'\t')
    f.write(str( E_cation_Rn   )+'\t')
    f.write(str( E_dication_Rd )+'\t')
    f.write(str( E_dication_Rc )+'\t')
    f.write(str( E_dication_Rn )+'\t')
    f.write(str( EPZ_neutral   )+'\t')
    f.write(str( EPZ_cation    )+'\t')
    f.write(str( EPZ_dication  ))
    f.write('\n')
    f.close()
    


def save_optimized_structure(n, method, species, label):
    """Write the converged geometry of a species. label: neutral, cation, dication."""
    write('structures/optimized/' + n + '_' + method + '_'+ label +'_bfgs.xyz', species)
        
    
# initialization
names = ['naphthalene', 'fluorene', 'anthracene', 'pyrene', 'tetracene', 'coronene', 'ovalene', 'dicoronylene']

methods = ['HF-3C', 'B97-3C', 'R2SCAN-3C', 'PBEH-3C', 'B3LYP-3C', 'WB97X-3C']

NPROCS = 32
MAXCORE = 4000                      # MB per process

# Molecules for which the Hueckel guess is necessary. Without it, the SCF of
# dicoronylene does not converge; on the smaller molecules the default guess of
# ORCA suffices.

NEEDS_HUECKEL = ['dicoronylene']


def orca_blocks(molecule, hftyp=None):
    """%scf block to pass to ORCA for this molecule.

    hftyp: 'ROHF' or 'RHF' to impose the type of wavefunction; None lets ORCA
    choose, that is RHF or UHF according to the multiplicity. The ROHF
    diagnostic of the HF-3C cations is done by a separate script
    (spin_contamination_HF-3C.py) and not here.
    """
    scf = []
    if molecule in NEEDS_HUECKEL:
        scf.append('guess hueckel')
    if hftyp is not None:
        scf.append(f'HFTyp {hftyp}')
    head = f'%scf {" ".join(scf)} end\n' if scf else ''
    return f'{head}%pal nprocs {NPROCS} end \n%maxcore {MAXCORE}'

max_force = 0.0005

fname = "results_"+time.strftime("%Y%m%d-%H%M%S")+".txt"

f = open(fname, "w")
f.close()

f = open(fname, "a")
f.write('Molecule'+'\t'+'Method'+'\t'+'E_neutral_Rn'+'\t'+'E_cation_Rc'+'\t'+'E_cation_Rn'+'\t'+'E_dication_Rd'+'\t'+'E_dication_Rc'+'\t'+'E_dication_Rn'+'\t'+'EPZ_neutral'+'\t'+'EPZ_cation'+'\t'+'EPZ_dication')
f.write('\n')
f.close()


# main loop
for method in methods:
    print()
    print("-----------------" + method + "-----------------")

    for n in names:
        print("             " + n)
        
        # Restart after an interruption: if the converged geometry already
        # exists, we resume from it and BFGS stops at the first step. Otherwise
        # we start from the initial structure of the astrochemical database.
        print("\nNeutral geometric optimization")        
        fname_opt = 'structures/optimized/' + n + '_' + method + '_neutral_bfgs.xyz'
        fname_ini = '../structures/initial/' + n + '.xyz'                              # from the astrochemistry database
        #print(fname_opt)
        #print(fname_ini)
        #print(os.getcwd())
        
        if os.path.isfile(fname_opt):
            print("Found neutral optimized structures:", fname_opt)
            neutral  = read(fname_opt)
        else:
            neutral  = read(fname_ini)  
            print("Found neutral initial structure: ", fname_ini)
        neutral.center()
        
        neutral.calc = ORCA(charge=0, mult=1,
                            label = n+'_neutral',
                            profile=profile,
                            # engrad is indispensable to the optimisations:
                            # ase.calculators.orca copies orcasimpleinput as it
                            # stands and never adds a gradient keyword. The
                            # forces are read back from orca.engrad, a file ORCA
                            # only writes if asked for it.
                            orcasimpleinput='engrad '+method,
                            orcablocks=orca_blocks(n),
                            )
        dyn = BFGS(neutral, trajectory='structures/optimized/' + n + '_' + method + '_neutral_bfgs.traj')
        dyn.run(fmax=max_force)
        save_optimized_structure(n, method, neutral, "neutral")
        
        E_neutral_Rn = neutral.get_potential_energy() #

        # Zero-point energy by finite differences of the forces: 6N + 1 gradient
        # evaluations. BEWARE: ase.vibrations caches each displacement in the
        # directory named below and, on a rerun, skips those that already exist,
        # WITHOUT checking that the geometry is the same. If a species is
        # reoptimised, its vib/ directory must be deleted, otherwise the ZPE
        # returned is that of the old geometry.
        print("Neutral vibrational calculation")
        vibr = 'vib/'+ 'vib_' + n + '_'+ method +'_neutral'   
        vib = Vibrations(neutral, name=vibr)
        vib.run()
        
        EPZ_neutral = vib.get_zero_point_energy()

        #----------- cation ---------
        print("\nCation geometric optimization")
        fname_opt = 'structures/optimized/' + n + '_' + method + '_cation_bfgs.xyz'
        fname_ini = '../structures/initial/' + n + '+.xyz'                              # from the astrochemistry database
        if os.path.isfile(fname_opt):
            print("Found cation optimized structures:", fname_opt)
            cation  = read(fname_opt)        
        else:
            #cation = neutral
            #print("Using neutral optimized structure")
            cation  = read(fname_ini)  
            print("Found cation initial structure: ", fname_ini)
        cation.center()
            
        #if method == "HF-3C":               # Keep for later, when the code is re-run for HF-3C with ROHF
        #    block = orca_block_ROHF         # At the moment, we optimize everything with RHF
        #else:
        #    block = orca_block
                
        cation.calc = ORCA(charge=1, mult=2,
                           label = n+'_cation',
                           profile=profile,
                           orcasimpleinput='engrad '+method,
                           orcablocks=orca_blocks(n),
                           )
        dyn = BFGS(cation, trajectory= 'structures/optimized/' + n + '_' + method + '_cation_bfgs.traj')
        dyn.run(fmax=max_force)
        save_optimized_structure(n, method, cation, "cation")
        
        E_cation_Rc = cation.get_potential_energy()

        print("cation vibrational calculation")
        vibr = 'vib/'+ 'vib_' + n  + '_'+ method +'_cation'
        vib = Vibrations(cation, name=vibr)
        vib.run()
        
        EPZ_cation = vib.get_zero_point_energy()

        # Cation at the geometry of the neutral. We reuse the neutral object,
        # whose positions are those converged above, changing only the charge
        # and the multiplicity of the calculator. No engrad here: a single point
        # suffices, the gradient would be paid for nothing.
        neutral.calc = ORCA(charge=1, mult=2,         #cation at the neutral geometry
                            label = n+'_cation_Rn',
                            profile=profile,
                            orcasimpleinput=method,
                            orcablocks=orca_blocks(n),
                            )
        E_cation_Rn = neutral.get_potential_energy()  # cation at the geometry of the neutral
        
        
        
        #block = orca_block


        #----------- dication -----------
        print("\nDication geometric optimization")
        # Multiplicity of the ground state of the dication, taken from the
        # theoretical spectral database of PAHs of Malloci et al. On this set,
        # only the coronene dication is a triplet; all the others are singlets.
        # The starting structures follow the same naming convention: '++t' for
        # the triplet, '++s' for the singlet.
        if n == 'coronene':
            fname_opt = 'structures/optimized/' + n + '_' + method + '_dication_bfgs.xyz'
            fname_ini = '../structures/initial/' + n + '++t.xyz'                              # from the astrochemistry database

        else:
            fname_opt = 'structures/optimized/' + n + '_' + method + '_dication_bfgs.xyz'
            fname_ini = '../structures/initial/' + n + '++s.xyz'                              # from the astrochemistry database
        
        if os.path.isfile(fname_opt):
            print("Found dication optimized structures:", fname_opt)
            dication  = read(fname_opt)
        else:
            #dication = cation
            #print("Using cation optimized structure")
            dication  = read(fname_ini)  
            print("Found dication initial structure: ", fname_ini)
        dication.center()

        if n == 'coronene':
            dication.calc = ORCA(charge=2, mult=3,
                                 label = n+'_dication',
                                 orcasimpleinput='engrad '+method,
                                 profile=profile,
                                 orcablocks=orca_blocks(n),
                                 )
        #elif n == 'coronene' and method == 'HF-3C':    #To used after RHF optimization
        #    dication.calc = ORCA(charge=2, mult=3,
        #                         label = n+'_dication',
        #                         orcasimpleinput=method,
        #                         profile=profile,
        #                         orcablocks=orca_block_ROHF,
        #                         )
        else:
            dication.calc = ORCA(charge=2, mult=1,
                                 label = n+'_dication',
                                 profile=profile,
                                 orcasimpleinput='engrad '+method,
                                 orcablocks=orca_blocks(n),
                                 )
        dyn = BFGS(dication, trajectory= 'structures/optimized/' + n + '_' + method + '_dication_bfgs.traj')
        dyn.run(fmax=max_force)
        save_optimized_structure(n, method, dication, "dication")
        
        E_dication_Rd = dication.get_potential_energy()
                
        print("dication vibrational calculation")
        vibr = 'vib/'+ 'vib_' + n + '_'+ method +'_dication'
        vib = Vibrations(dication, name=vibr)
        vib.run()
        
        EPZ_dication = vib.get_zero_point_energy()
    

        if n == 'coronene':
            cation.calc = ORCA(charge=2, mult=3,        #dication at the cation geometry
                               label = n+'_dication_Rc',
                               profile=profile,
                               orcasimpleinput=method,
                               orcablocks=orca_blocks(n),
                               )
        #elif n == 'coronene' and method == 'HF-3C':
        #   cation.calc = ORCA(charge=2, mult=3,        #dication at the cation geometry
        #                      label = n+'_dication_Rc',
        #                      profile=profile,
        #                      orcasimpleinput=method,
        #                      orcablocks=orca_block_ROHF,
        #                      )
        else:
            cation.calc = ORCA(charge=2, mult=1,       #dication at the cation geometry
                               label = n+'_dication_Rc',
                               profile=profile,
                               orcasimpleinput=method,
                               orcablocks=orca_blocks(n),
                               )
        E_dication_Rc = cation.get_potential_energy()   # dication at the geometry of the cation
        

        if n == 'coronene':
            neutral.calc = ORCA(charge=2, mult=3,       #dication at the neutral geometry
                                label = n+'_dication_Rn',
                                profile=profile,
                                orcasimpleinput=method,
                                orcablocks=orca_blocks(n),
                                )
        #elif n == 'coronene' and method == 'HF-3C'
        #    neutral.calc = ORCA(charge=2, mult=3,       #dication at the neutral geometry
        #                        label = n+'_dication_Rn',
        #                        profile=profile,
        #                        orcasimpleinput=method,
        #                        orcablocks=orca_block_ROHF,
        #                        )
        else:
            neutral.calc = ORCA(charge=2, mult=1,       #dication at the neutral geometry
                                label = n+'_dication_Rn',
                                profile=profile,
                                orcasimpleinput=method,
                                orcablocks=orca_blocks(n),
                                )
        E_dication_Rn = neutral.get_potential_energy()   # dication at the geometry of the neutral
        

        print(E_neutral_Rn, E_cation_Rc, E_cation_Rn, E_dication_Rd, E_dication_Rc, E_dication_Rn, EPZ_neutral, EPZ_cation, EPZ_dication )

        saveenergies(n, method, fname, E_neutral_Rn, E_cation_Rc, E_cation_Rn, E_dication_Rd, E_dication_Rc, E_dication_Rn, EPZ_neutral, EPZ_cation, EPZ_dication)

        print("\n\n")

        os.system('rm orca.gbw')  # clean up command to avoid discrepancy between input geometry and actual geometry
        
    print("End of method ", method, "\n")

    
print("\n \n -------- END OF CALCULATIONS-------- \n")

# Ionization energies of polycyclic aromatic hydrocarbons — benchmark dataset

Supporting dataset for *Ionization Energies of Polycyclic Aromatic Hydrocarbons:
Benchmarking Semi-Empirical, Tight Binding and Composite Methods*, A. Giuliani and
H. R. Hrodmarsson, ChemPhysChem.

This repository contains the optimized geometries, the raw total energies and the analysis scripts from which every number of the article and of its Supporting Information can be regenerated.

## Contents

```
geometries/   312 optimized structures, XYZ format, Angstrom
data/         raw total energies and zero-point energies, eV
scripts/      generation and analysis scripts
```

### geometries/

Naming: `<molecule>_<method>_<species>_bfgs.xyz`

- **molecules** (8): naphthalene, fluorene, anthracene, pyrene, tetracene,  coronene, ovalene, dicoronylene
- **methods** (13): DFTB2, DFTB3, IPEA1, PM6, PM7, RM1, g-xTB, HF-3C, B97-3C, R2SCAN-3C, PBEH-3C, B3LYP-3C, WB97X-3C
- **species** (3): neutral, cation, dication

8 × 13 × 3 = 312 files. Coordinates are in Angstrom. All structures were relaxed with the ASE BFGS optimizer to a maximum force of 5 × 10⁻⁴ eV/Å.

Charges and multiplicities are implicit in the species label: neutral = 0/singlet, cation = +1/doublet, dication = +2/singlet, except for the coronene dication whose ground state is a triplet.

### data/

Tab-separated, one row per molecule and method. Energies in eV.

| file                           | content                                                                                  |
|--------------------------------|------------------------------------------------------------------------------------------|
| `results_semi-empirical.txt`   | DFTB2, DFTB3, IPEA1, PM6, PM7, RM1                                                       |
| `results_g-xTB.txt`            | g-xTB                                                                                    |
| `results_composite.txt`        | HF-3C, B97-3C, R2SCAN-3C, PBEH-3C, B3LYP-3C, WB97X-3C                                    |
| `spin_contamination_HF-3C.txt` | ⟨S²⟩ and ROHF/UHF energy differences for the HF-3C open-shell species                    |
| `timings_coronene.txt`         | computational cost of one energy and one energy-plus-gradient evaluation, all 13 methods |

Columns of the three `results_*.txt` files:

```
Molecule  Method
E_neutral_Rn    total energy of the neutral   at the neutral geometry
E_cation_Rc     total energy of the cation    at the cation geometry
E_cation_Rn     total energy of the cation    at the neutral geometry
E_dication_Rd   total energy of the dication  at the dication geometry
E_dication_Rc   total energy of the dication  at the cation geometry
E_dication_Rn   total energy of the dication  at the neutral geometry
EPZ_neutral     zero-point energy of the neutral
EPZ_cation      zero-point energy of the cation
EPZ_dication    zero-point energy of the dication
```

Adiabatic and vertical ionization energies follow from these six energies, with the zero-point corrections added where appropriate:

```
aIE1 = E_cation_Rc  - E_neutral_Rn        vIE1 = E_cation_Rn  - E_neutral_Rn
aIE2 = E_dication_Rd - E_cation_Rc        vIE2 = E_dication_Rc - E_cation_Rc
aIEd = E_dication_Rd - E_neutral_Rn       vIEd = E_dication_Rn - E_neutral_Rn
```

Two conventions used throughout the article are applied by the analysis script and must be reproduced to obtain the published numbers:

1. The zero-point energies of g-xTB are absent, as vibrational frequencies could not be computed with the preliminary version used. All g-xTB ionization energies are therefore reported without zero-point correction.
2. xTB-IPEA requires the empirical correction of 4.8455 eV per ionization step prescribed by the authors of the method; it is subtracted once for the first and second ionization energies and twice for the double ionization energy.

### scripts/

| file                          | purpose                                                                                             |
|-------------------------------|-----------------------------------------------------------------------------------------------------|
| `benchmark-approximate.py`    | production script for the semi-empirical and tight-binding methods (ASE + DFTB+ + MOPAC)            |
| `benchmark-gxtb.py`           | production script for g-xTB (ASE + the preliminary g-xTB binary)                                    |
| `benchmark-3C.py`             | production script for the composite methods (ASE + ORCA)                                            |
| `spin_contamination_HF-3C.py` | ⟨S²⟩ and ROHF/UHF diagnostic for the HF-3C open-shell species                                       |
| `check_results.py`            | checks the physical consistency of the deposited energies                                           |
| `analysis.jl`                 | regenerates Tables 1 to 3, S2, S3, S4, S6 and S7 and both figures from `data/`, and writes them out |
| `timing_benchmark.py`         | measures the computational cost of each method; this is what produced `data/timings_coronene.txt`   |
| `gxtb.py`                     | ASE calculator for g-xTB, required by `benchmark-gxtb.py` and `timing_benchmark.py`                 |

Three traps were met while timing the methods, each of which silently produced a plausible but wrong result, and each of which is documented in the comments of `timing_benchmark.py`: ASE never adds a gradient keyword to the ORCA input, so the
forces must be requested explicitly; the ASE MOPAC and g-xTB calculators request the gradient by default, so the energy-only step has to disable it; and `Calculator.reset()` does not exist on the calculator family ORCA belongs to. Anyone reusing the script should keep those in mind.

`timing_benchmark.py` measures one energy and one energy-plus-gradient evaluation, always on a single core, and derives the cost of a Hessian from the latter as (6N + 1) × t_grad rather than measuring it. Earlier versions could also run in parallel and could time the frequencies directly; both were removed, the first because wall-clock times obtained on different numbers of cores do not compare, the second because measuring the Hessian only repeats the gradient for several hours per method. A file of measured frequencies for naphthalene, obtained with that earlier version, was deposited with a previous release and has been withdrawn, the present script no longer being able to produce it.

The comments of the scripts are in English, their runtime messages in French.

### Running the production scripts

The three `benchmark-*.py` scripts are run from a directory holding `structures/initial/`, the starting geometries taken from the astrochemical database, named `<molecule>.xyz`, `<molecule>+.xyz` and `<molecule>++s.xyz`. They write the relaxed geometries to `structures/optimized/`, the finite-difference gradient caches to `vib/`, and one tab-separated `results_<timestamp>.txt`.

Code locations default to `~/dftb+`, `~/g-xtb` and `~/local`, expanded at run time because neither DFTB+ nor g-xTB interprets a tilde. Each is overridable, so an installation elsewhere needs no edit to the scripts:

```
DFTB_SK_HOME    Slater-Koster sets, expected to hold mio-1-1/ and 3ob-3-1/
GXTBHOME        g-xTB parameter files (.gxtb, .basisq, .eeq)
GXTB_COMMAND    g-xTB executable
ORCA_BIN        ORCA executable, full path
LOCAL_CODES     root used by timing_benchmark.py; see its header for the rest
```

`ORCA_BIN` has no usable default and must be set: ORCA is rarely on the PATH, and on most distributions `/usr/bin/orca` is the GNOME screen reader rather than the quantum chemistry program.

Two prerequisites are not satisfied by a stock installation. `gxtb.py` must be installed as `ase/calculators/gxtb.py` for `benchmark-gxtb.py` to import it. And
`benchmark-approximate.py` passes `mult='TRIPLET UHF'` to the ASE MOPAC calculator, a keyword the distributed ASE does not have: it needs a local addition to `ase/calculators/mopac.py` providing the `ms` and `mult` parameters and forcing the `DISP` keyword so that the total energy is reported. Without it the coronene dication would be computed closed-shell and the keyword silently ignored.

### Reproducing the tables and the figures

`analysis.jl` regenerates everything the article derives from the raw energies, and writes it to `output/`:

```
julia scripts/analysis.jl
```

| written to `output/`                                 | content                                                               |
|------------------------------------------------------|-----------------------------------------------------------------------|
| `table1_IE1.tex`, `table2_IE2.tex`, `table3_IEd.tex` | first, second and double adiabatic ionization energies, ZPE-corrected |
| `tableS2_ionization.tex`                             | adiabatic energies with and without ZPE, and vertical energies        |
| `tableS3_metrics.tex`                                | MSE, MAE, RMSE and σ for the thirteen methods                         |
| `tableS4_zpe.tex`                                    | zero-point energies of the neutral, and ΔZPE                          |
| `tableS6_topology.tex`                               | deviations along the linear-acene and pericondensed subseries         |
| `tableS7_cost.tex`                                   | cost of one energy evaluation, relative to DFTB2                      |
| `box.png`                                            | Figure 2, box plots of the errors                                     |
| `gToC.png`                                           | the graphical table of contents                                       |

The `.tex` files hold the rows of the present work only, and the numbers alone: the tabular preamble, the thin space of the thousands and the `$-$` of the negative values belong to the manuscript. The literature rows of Tables 1 to 3 are quoted from the cited references and are not reproduced here, no more than the xTB-GFN2 row of Table S4, which is read from the COMPAS database. Both figures come out byte for byte identical to the published ones.

Tables S1 and S5 are not written. S1 is the content of `data/results_*.txt` as it stands, and S5 that of `data/spin_contamination_HF-3C.txt`, the two differences it reports having been taken by the script that produced that file. Neither involves an operation on the deposited numbers.

One caveat on Table S4. ΔZPE is computed here at full precision and rounded once, at printing. The published table does not do so throughout: its DFTB2, xTB-IPEA, PM6, PM7 and RM1 rows are the difference of the two zero-point energies already rounded to three decimals, which shifts fifteen of its ninety-six values by 0.001 eV. The DFTB3 row and the six composite rows are at full precision, as here.

It needs Julia with:

```
julia -e 'using Pkg; Pkg.add(["DataFrames","CSV","StatsPlots","Plots"])'
```

### Checking the raw energies

`check_results.py` verifies the physical consistency of the deposited energies. It takes no argument and checks the three `results_*.txt` files; to check a fresh campaign before deposit, add its path to the `DEPOSITED` list at the top of the script:

```
python3 scripts/check_results.py
```

The three relations it tests are variational bounds — a relaxed geometry minimises the energy of its own species, so `E_cation(R_c)` and `E_dication(R_d)` must lie below the same species at any other geometry:

```
E_cation(R_c)   ≤ E_cation(R_n)
E_dication(R_d) ≤ E_dication(R_c)
E_dication(R_d) ≤ E_dication(R_n)
```

Violating one means an optimisation that has not converged, a charge or a spin state wrongly set, or a geometry read from the wrong file. The script returns 1 in that case, 0 otherwise, so that it can guard a deposit. On the 104 deposited rows it returns 0.

Only what would denounce a defect is tested. Orderings that depend on the relative position of the three geometries, without the standing of a theorem, are not: `E_dication(R_c) ≤ E_dication(R_n)` holds on all but two rows of the deposit — coronene in B97-3C and in R2SCAN-3C, by 2 and 5 meV — coronene being the one molecule of the set whose dication is a triplet, and whose doublet cation is Jahn-Teller distorted in a direction that does not lead to it. Those are not errors, and are not flagged.

The ionization energies and the statistical metrics are not recomputed here; `analysis.jl` owns them.

## Software

- ASE 3.23
- ORCA 6.1.0 (composite methods)
- DFTB+ 23.1 (DFTB2, DFTB3, xTB-IPEA)
- OpenMOPAC 23.0.1 (PM6, PM7, RM1)
- g-xTB, preliminary g-xTB(EEQ BC) build of 12 August 2025
 
The g-xTB executable is a preliminary version that is available here: https://github.com/grimme-lab/g-xtb .

## Citing

Please cite the article above. If you use the dataset itself, a reference to this repository is appreciated.

## Licence

Copyright © 2026 A. Giuliani. Two licences apply, one to the content and one to
the code, because CC BY is written for data and documents and says nothing useful
about patents or object code.

**`geometries/` and `data/`** — Creative Commons Attribution 4.0 International
(CC BY 4.0). Full legal code in `LICENSE`, summary at
https://creativecommons.org/licenses/by/4.0/ . You are free to share and adapt
them, for any purpose including commercially, provided you give appropriate
credit.

**`scripts/`** — MIT, in `scripts/LICENSE`, with one exception. `gxtb.py` was
written using the `dftb.py` calculator of ASE as a template, and ASE is
distributed under LGPL-2.1-or-later; that file is therefore covered by
LGPL-2.1-or-later rather than by MIT. The other seven scripts are original work.

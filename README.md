# Ionization energies of polycyclic aromatic hydrocarbons — benchmark dataset

Supporting dataset for *Ionization Energies of Polycyclic Aromatic Hydrocarbons: Benchmarking Semi-Empirical, Tight Binding and Composite Methods, A. Giuliani and H. R. Hrodmarsson, ChemPhysChem 2026, 27, e70622. https://doi.org/10.1002/cphc.70622.

This repository contains the optimized geometries, the raw total energies and the analysis scripts from which every number of the article and of its Supporting Information can be regenerated.

## Contents

```
geometries/   313 optimized structures, XYZ format, Angstrom
data/         raw total energies and zero-point energies, eV
scripts/      generation and analysis scripts
```

### geometries/

Naming: `<molecule>_<method>_<species>_bfgs.xyz`

- **molecules** (8): naphthalene, fluorene, anthracene, pyrene, tetracene,  coronene, ovalene, dicoronylene
- **methods** (13): DFTB2, DFTB3, IPEA1, PM6, PM7, RM1, g-xTB, HF-3C, B97-3C, R2SCAN-3C, PBEH-3C, B3LYP-3C, WB97X-3C
- **species** (3): neutral, cation, dication

8 × 13 × 3 = 312 files, plus the one structure described just below. Coordinates are in Angstrom. All structures were relaxed with the ASE BFGS optimizer to a maximum force of 5 × 10⁻⁴ eV/Å.

Charges and multiplicities are implicit in the species label: neutral = 0/singlet, cation = +1/doublet, dication = +2/singlet, except for the coronene dication whose ground state is a triplet.

One further structure sits outside that grid:

```
circumcircumcoronene_DFTB3_neutral_bfgs.xyz    C96H24, 120 atoms
```

It belongs to the cost measurement only, not to the benchmark set: no ionization energy is reported for it, and only its neutral is needed, since the cost of an evaluation does not depend on the charge. It is a zigzag hexagonal flake of the coronene family (37 hexagons, honeycomb construction, C–C 1.42 Å then H at 1.09 Å), relaxed with DFTB3/3ob in 24 BFGS steps to the same 5 × 10⁻⁴ eV/Å: planar, D6h, C–C 1.364 to 1.444 Å. The same construction gives coronene back for n = 2. All thirteen methods were timed from this one geometry.

### data/

Tab-separated, one row per molecule and method. Energies in eV.

| file                           | content                                                                                  |
|--------------------------------|------------------------------------------------------------------------------------------|
| `results_semi-empirical.txt`   | DFTB2, DFTB3, IPEA1, PM6, PM7, RM1                                                       |
| `results_g-xTB.txt`            | g-xTB                                                                                    |
| `results_composite.txt`        | HF-3C, B97-3C, R2SCAN-3C, PBEH-3C, B3LYP-3C, WB97X-3C                                    |
| `spin_contamination_HF-3C.txt` | ⟨S²⟩ and ROHF/UHF energy differences for the HF-3C open-shell species                    |
| `timings_coronene.txt`         | cost of one energy and one energy-plus-gradient evaluation on coronene, all 13 methods   |
| `timings_size_scaling.txt`     | cost of one energy evaluation at four sizes, 36 to 120 atoms, all 13 methods             |

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

The statistics rest on the molecules for which an experimental value is retained: seven for the first ionization energy (dicoronylene has none), six for the second (fluorene has none, and the naphthalene value, derived from charge-stripping data, is set aside) and six for the double ionization energy (ovalene and dicoronylene have none). The submitted version used five for the double ionization energy, anthracene having been left out by mistake; the revised version, and the present scripts, include it. The experimental values and their quoted uncertainties are the arrays `IE*_EXP` and `IE*_UNC` at the top of `analysis.jl`.

#### The two cost files

They were measured a month and a half apart, under conditions that are not the same, and they are not interchangeable.

`timings_coronene.txt` (29 July 2026) was the source of Table S7 of the submitted version, which the revised version replaces with Figure S1. It is kept because it is the only measurement of the gradient, on which the discussion of the cost of the complete workflow rests. It holds `t_energy_s` and `t_energy_grad_s` on coronene alone, and the Hessian cost derived from the latter as (6N + 1) × t_grad. Its single points are **warm**: ORCA re-reads its own `orca.gbw` between the trials of `--repeat` unless told not to, so that from the second trial on the SCF starts from a density already converged for the same geometry. `t_energy_grad_s`, the quantity the file is for, is little affected, the gradient being a fixed post-SCF cost that a good initial density does not speed up.

`timings_size_scaling.txt` (15–16 September 2026) is the source of Figure S1, the cost-against-size figure. Its single points are **cold**: `NoAutoStart` is now added to every ORCA input, so each trial starts from the model density, as any real workflow does. It holds `t_sp` alone — measured with `--no-grad`, no gradient and no Hessian — at four sizes:

```
coronene              C24H12    36 atoms
ovalene               C32H14    46 atoms
dicoronylene          C48H20    68 atoms
circumcircumcoronene  C96H24   120 atoms
```

together with `alpha`, the slope of ln t against ln N over the four sizes, `alpha_corr`, the same after subtraction of the fixed start-up overhead of 0.034 s that only matters for the fast methods, and `ratio_120_36`. The three derived columns follow from the four times and can be recomputed from them.

The coronene column was remeasured rather than taken from July, for two reasons: it is cold where July was warm, and ORCA itself ran about twice as slowly in July as in September for an unidentified reason, while the fast methods were unchanged (DFTB2 0.084 s against 0.082 s). A single machine state therefore carries the whole table. The two effects do not cancel, and the July and September coronene columns differ by factors of 0.8 to 2.4 depending on the method.

The last sizes of each stream were measured with less of the machine busy, the HF-3C control runs putting the cost of concurrency at up to 25 %; that lowers the exponents by at most 0.2.

### scripts/

| file                          | purpose                                                                                             |
|-------------------------------|-----------------------------------------------------------------------------------------------------|
| `benchmark-approximate.py`    | production script for the semi-empirical and tight-binding methods (ASE + DFTB+ + MOPAC)            |
| `benchmark-gxtb.py`           | production script for g-xTB (ASE + the preliminary g-xTB binary)                                    |
| `benchmark-3C.py`             | production script for the composite methods (ASE + ORCA)                                            |
| `spin_contamination_HF-3C.py` | ⟨S²⟩ and ROHF/UHF diagnostic for the HF-3C open-shell species                                       |
| `check_results.py`            | checks the physical consistency of the deposited energies                                           |
| `analysis.jl`                 | regenerates Tables 1 to 3, S2, S3, S4 and S6, Figure 2 and the graphical ToC from `data/`          |
| `resampling.jl`               | resolution of the statistical metrics: regenerates Tables S7 and S8 (sample size, experimental uncertainty, combined) |
| `timing_benchmark.py`         | measures the computational cost of each method; produced both `data/timings_*.txt`                  |
| `figure_cost_vs_size.jl`      | draws Figure S1, the cost-against-size figure, from `data/timings_size_scaling.txt`                 |
| `gxtb.py`                     | ASE calculator for g-xTB, required by `benchmark-gxtb.py` and `timing_benchmark.py`                 |

Three traps were met while timing the methods, each of which silently produced a plausible but wrong result, and each of which is documented in the comments of `timing_benchmark.py`: ASE never adds a gradient keyword to the ORCA input, so the
forces must be requested explicitly; the ASE MOPAC and g-xTB calculators request the gradient by default, so the energy-only step has to disable it; and `Calculator.reset()` does not exist on the calculator family ORCA belongs to. Anyone reusing the script should keep those in mind.

`timing_benchmark.py` measures one energy and one energy-plus-gradient evaluation, always on a single core, and derives the cost of a Hessian from the latter as (6N + 1) × t_grad rather than measuring it. Earlier versions could also run in parallel and could time the frequencies directly; both were removed, the first because wall-clock times obtained on different numbers of cores do not compare, the second because measuring the Hessian only repeats the gradient for several hours per method. A file of measured frequencies for naphthalene, obtained with that earlier version, was deposited with a previous release and has been withdrawn, the present script no longer being able to produce it.

Three more things were learned in September 2026, when the measurement was extended to four sizes, and each is documented in the comments of the scripts:

- ORCA restarts from its own converged density unless `NoAutoStart` is given, which makes every trial after the first a warm one. The keyword is now always added, and `--no-grad` was added with it, to measure `t_sp` alone at sizes where a gradient is not wanted.
- With ORCA 6.1.1 the default model-potential guess aborts for ωB97X-3C on dicoronylene inside `orca_guess` (*error with the V-matrix from which there is no recovery*) while it goes through on coronene. The extended Hückel guess goes through, and is imposed on that method at every size so that the times stay comparable. The other composite methods keep the ORCA default.
- The g-xTB binary stops printing its per-atom table beyond about a hundred atoms, on which `gxtb.py` raised `UnboundLocalError` although the calculation had run and the energy was correct. `read_charges_dipole` now returns `None` when the table is absent.

Concurrent streams measuring the same method on the same molecule used to write the same working directory and make each other fail; `--tag` now labels both the working directory and the output file.

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
| `box.png`                                            | Figure 2, box plots of the errors                                     |
| `gToC.png`                                           | the graphical table of contents                                       |

The `.tex` files hold the rows of the present work only, and the numbers alone: the tabular preamble, the thin space of the thousands and the `$-$` of the negative values belong to the manuscript. The literature rows of Tables 1 to 3 are quoted from the cited references and are not reproduced here, no more than the xTB-GFN2 row of Table S4, which is read from the COMPAS database. Both figures come out byte for byte identical to the published ones.

The cost-against-size figure is drawn by a script of its own, since it rests on `data/timings_size_scaling.txt` and on nothing `analysis.jl` reads:

```
julia scripts/figure_cost_vs_size.jl
```

| written to `output/`                                     | content                                                        |
|----------------------------------------------------------|----------------------------------------------------------------|
| `figure_cost_vs_size.png`, `figure_cost_vs_size.pdf`     | Figure S1, single-point wall time against number of atoms, log-log, 300 dpi |

It needs `Plots` alone, already installed by the line above. Circles mark the fast methods and squares the composite ones, the colours are those of the other figures, and the number in brackets after each label is the effective exponent α, the raw slope over the four sizes (the `alpha` column of the data file). It too comes out byte for byte identical to the published figure.

### Resolution of the statistical metrics

The metrics rest on six or seven molecules, and the reference values carry experimental uncertainties of up to 0.2 eV. `resampling.jl` estimates how precisely the mean absolute error (MAE) of each method is known, as 95 % half-widths (eV):

```
julia scripts/resampling.jl
```

- **(A) experimental uncertainty alone**: each reference value is replaced by a draw from N(e_i, u_i²), the calculated values being held fixed, and the MAE is recomputed;
- **(B) sample size alone**: the molecules are drawn with replacement among themselves (bootstrap) and the MAE is recomputed;
- **(C) both combined**: the two sources being independent, C = √(A² + B²).

(A) and (B) use 2 × 10⁵ draws each, with a fixed seed, so the output is reproducible; changing the seed moves the values by at most 0.01 eV. The script includes `analysis.jl` for the data, the ionization energies and the uncertainties, so nothing is duplicated.

| written to `output/`         | content                                                                                   |
|------------------------------|-------------------------------------------------------------------------------------------|
| `tableS7_resolution.tex`     | Table S7, the three half-widths, median over the thirteen methods                        |
| `tableS8_combined.tex`       | Table S8, the MAE of each method and its three half-widths                               |
| `resolution_full.txt`        | the same per method, with the RMSE and the analytic standard error s/√n as a check on (B) |

With the median combined half-widths U (0.05, 0.11 and 0.25 eV), two methods can be told apart only if their MAE differ by more than about √2 U, that is 0.1, 0.15 and 0.3 eV for the first, second and double ionization energies: the thresholds used in the article.

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

The energies were produced with those versions. The cost measurements were made on one machine, a Xeon E5-1650 v2 under Ubuntu 24.04, with ASE 3.28 and ORCA 6.1.1 in September 2026; the absolute times are proper to it, the ratios between methods much less so.
 
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

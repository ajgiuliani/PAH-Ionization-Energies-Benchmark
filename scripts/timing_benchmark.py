#!/usr/bin/env python3
"""
Computational cost of the methods of the PAH benchmark.

The production logs do not allow a publishable table of costs to be drawn up:
they cover sessions that were interrupted and resumed, so that the intervals
between blocks include dead time (up to 22 h on some blocks). This script redoes
the measurement under controlled conditions.

PURPOSE
-------
To give an ORDER OF MAGNITUDE of the cost of each method.

Protocol: one single molecule, coronene (C24H12, 36 atoms), the median size of
the benchmark set. For each method, two quantities are timed separately

  1. t_sp   : one energy evaluation alone;
  2. t_grad : one energy and gradient evaluation.


The cost of a Hessian follows from t_grad: ase.vibrations proceeds by finite
differences of the forces, that is 6 N_atoms + 1 gradient evaluations. It is
therefore estimated and printed as a derived quantity, never measured: measuring
it would only repeat t_grad, for several hours per method on coronene.

One option measures what is otherwise estimated:

  --with-opt    measures a complete optimisation; its time depends on the
                starting point and on the optimiser, and is to be read as an
                indication.

Two options were added for the size-scaling campaign of 15-16 September 2026:

  --no-grad     measures t_sp alone; t_grad and the estimate of the Hessian are
                then nan. This is what gives a size column (ovalene,
                dicoronylene, circumcircumcoronene) to a table that compares
                single points only.
  --tag         a label inserted in the name of the output file and of the
                working directory, so that two concurrent streams measuring the
                same method on the same molecule cannot write the same files.

COLD AND WARM SINGLE POINTS
---------------------------
ORCA re-reads its own orca.gbw (AutoStart, MOREAD) unless told not to, so that
from the second trial of --repeat on, the SCF starts from a density already
converged for the SAME geometry, which no real workflow ever does. NoAutoStart
is now added to every ORCA input, and every trial therefore starts from the
model density: that is the cost of a COLD single point. The coronene campaign of
29 July 2026, which produced data/timings_coronene.txt, was measured WITHOUT the
keyword, that is warm; data/timings_size_scaling.txt is cold. The two files are
not interchangeable, and the difference reaches a factor of 4 on B97-3C.

The settings of the calculators reproduce exactly those of the production
scripts (benchmark-approximate.py for the semi-empirical and tight-binding
methods, benchmark-gxtb.py for g-xTB, benchmark-3C.py for the composite ones).

RESOURCES AND COMPARABILITY
---------------------------
The methods do not parallelise in the same way. Comparing wall-clock times
obtained on different numbers of cores would therefore make no sense, and the
number of cores is not a parameter of the protocol: everything is forced onto a
single core. That is what gives a ratio of intrinsic cost between methods, and
that is how the published table was obtained.

Usage:
    python3 timing_benchmark.py                  # all the methods
    python3 timing_benchmark.py DFTB3 PM7        # a subset
    python3 timing_benchmark.py --with-opt       # with the optimisation
    python3 timing_benchmark.py --molecule pyrene
    python3 timing_benchmark.py --no-grad --molecule dicoronylene
"""

import argparse
import os
import shutil
import sys
import time

from ase.io import read
from ase.optimize import BFGS

# ---------------------------------------------------------------- configuration

MAX_FORCE = 0.0005               # same criterion as the production campaign


LOCAL = os.path.expanduser(os.environ.get("LOCAL_CODES", "~/local"))

ORCA_BIN = os.environ.get("ORCA_BIN", f"{LOCAL}/orca/orca")
DFTB_PREFIX_MIO = os.environ.get(
    "DFTB_PREFIX_MIO", f"{LOCAL}/SlaterKosterFiles/mio-1-1/")
DFTB_PREFIX_3OB = os.environ.get(
    "DFTB_PREFIX_3OB", f"{LOCAL}/SlaterKosterFiles/3ob-3-1/")
GXTB_BIN = os.environ.get("GXTB_BIN", f"{LOCAL}/g-xtb/binary/gxtb")

# gxtb only computes the gradient if it is passed -grad, and that gradient
# costs about 140 times the energy alone on coronene (0.4 s against 57 s). But
# the gxtb.py calculator adds -grad in ALL its default branches: without an
# explicit command=, the "energy alone" step was therefore computing the
# gradient too, and t_sp was a second t_grad. Both commands are thus supplied
# and make_calc chooses according to the step. gxtb.py is not modified: its
# __init__ keeps each default behind an "if command is None".
GXTB_COMMAND = os.environ.get(
    "GXTB_COMMAND", f"{GXTB_BIN} -c coord -grad > PREFIX.out")
GXTB_COMMAND_ENERGY = os.environ.get(
    "GXTB_COMMAND_ENERGY", f"{GXTB_BIN} -c coord > PREFIX.out")

# gxtb reads its parameters (.gxtb, .eeq, .basisq) through GXTBHOME; the binary
# does not start if that variable is unset, hence the default below.
os.environ.setdefault("GXTBHOME", f"{LOCAL}/g-xtb/parameters")

# The ASE calculator for g-xTB does not belong to ASE: it is a local module
# (gxtb.py) that must accompany this script on the compute machine. It is
# imported from ASE if it has been installed there, otherwise from the path
# given by GXTB_ASE_MODULE, otherwise from next to this script, otherwise from
# the neighbouring g-xTB directory.
_HERE = os.path.dirname(os.path.abspath(__file__))
GXTB_ASE_MODULE = os.environ.get(
    "GXTB_ASE_MODULE",
    next((p for p in (os.path.join(_HERE, "gxtb.py"),
                      os.path.join(_HERE, "g-xTB", "gxtb.py"))
          if os.path.isfile(p)), os.path.join(_HERE, "gxtb.py")))


def load_gxtb_class():
    try:
        from ase.calculators.gxtb import Gxtb
        return Gxtb
    except ImportError:
        pass
    if os.path.isfile(GXTB_ASE_MODULE):
        import importlib.util
        spec = importlib.util.spec_from_file_location("gxtb_local", GXTB_ASE_MODULE)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod.Gxtb
    return None

SEMI_EMPIRICAL = ["PM6", "PM7", "RM1"]
TIGHT_BINDING = ["DFTB2", "DFTB3", "IPEA1", "g-xTB"]
COMPOSITE = ["HF-3C", "B97-3C", "R2SCAN-3C", "PBEH-3C", "B3LYP-3C", "WB97X-3C"]
ALL_METHODS = TIGHT_BINDING[:3] + SEMI_EMPIRICAL + ["g-xTB"] + COMPOSITE


# ----------------------------------------------------------------- calculators

def _dftb_common(charge):
    """Options common to the three DFTB+ Hamiltonians, taken from production.

    UseOmpThreads is No: the measurement is made on a single core, whereas
    production opens it.
    """
    return dict(
        Hamiltonian_Charge=f"{charge}.0",
        Hamiltonian_MaxAngularMomentum_="",
        Hamiltonian_MaxAngularMomentum_C='"p"',
        Hamiltonian_MaxAngularMomentum_H='"s"',
        ParserOptions_="",
        ParserOptions_ParserVersion="10",
        Parallel="{ UseOmpThreads = No }",
    )


def make_calc(method, workdir, charge=0, forces=False):
    """ASE calculator configured as in the production scripts.

    forces: ask for the gradient on top of the energy. Concerns ORCA only.
    ase.calculators.orca copies orcasimpleinput as it stands into the .inp and
    never adds a gradient keyword, whatever property is requested
    (ase/io/orca.py: write_orca); the forces, for their part, are read back from
    orca.engrad, a file ORCA only writes if asked for it. Without EnGrad,
    get_forces() raises PropertyNotImplementedError AFTER having run ORCA:
    t_grad would then equal t_sp and the estimate of the Hessian would follow
    from it, plausible and wrong. The keyword is only added for the steps that
    need the gradient, so that t_sp remains an energy alone.
    """
    os.makedirs(workdir, exist_ok=True)

    if method in COMPOSITE:
        from ase.calculators.orca import ORCA, OrcaProfile
        # NoAutoStart: without it ORCA re-reads its own orca.gbw (AutoStart,
        # MOREAD) from the second trial of --repeat on, and the SCF converges in
        # two cycles on a density already converged for the SAME geometry, which
        # no real workflow does. Every trial therefore starts from the model
        # density: this is the cost of a cold single point. The coronene campaign
        # of 29 July 2026 was measured WITHOUT this keyword, that is warm.
        #
        # HUECKEL for WB97X-3C: with ORCA 6.1.1 the model-potential guess (the
        # default) aborts on dicoronylene inside orca_guess ("error with the
        # V-matrix from which there is no recovery", shark_setup_ri.cpp;
        # reproduced in serial and on two MPI ranks, 15 September 2026) while it
        # goes through on coronene. The extended Hueckel guess goes through; it
        # is imposed at ALL sizes so that the times remain comparable. The other
        # composite methods keep the ORCA default (B97-3C already uses EHT).
        guess = " HUECKEL" if method == "WB97X-3C" else ""
        return ORCA(charge=charge, mult=1,
                    directory=workdir,
                    profile=OrcaProfile(command=ORCA_BIN),
                    orcasimpleinput=(f"{method} NoAutoStart{guess} EnGrad" if forces
                                     else f"{method} NoAutoStart{guess}"),
                    orcablocks="%pal nprocs 1 end \n%maxcore 4000")

    if method in SEMI_EMPIRICAL:                      # PM6, PM7, RM1
        # The same trap as for g-xTB: ase.calculators.mopac fixes
        # task='1SCF GRADIENTS' in default_parameters, so that the two steps
        # were writing an identical .mop and measuring the same calculation
        # twice. The gradient only costs 0.016 s here out of 0.075 s of SCF (on
        # one thread), but that is enough for t_grad to come out sometimes
        # BELOW t_sp, which is impossible for a real gradient.
        from ase.calculators.mopac import MOPAC
        return MOPAC(method=method, charge=float(charge),
                     task="1SCF GRADIENTS" if forces else "1SCF",
                     label=os.path.join(workdir, "mopac"))

    if method == "DFTB2":
        from ase.calculators.dftb import Dftb
        return Dftb(label=os.path.join(workdir, "dftb"),
                    Hamiltonian_="DFTB",
                    Hamiltonian_SCC="Yes",
                    Hamiltonian_SlaterKosterFiles_="Type2FileNames",
                    Hamiltonian_SlaterKosterFiles_Prefix=DFTB_PREFIX_MIO,
                    Hamiltonian_SlaterKosterFiles_Separator='"-"',
                    Hamiltonian_SlaterKosterFiles_Suffix='".skf"',
                    **_dftb_common(charge))

    if method == "DFTB3":
        from ase.calculators.dftb import Dftb
        return Dftb(label=os.path.join(workdir, "dftb"),
                    Hamiltonian_="DFTB",
                    Hamiltonian_SCC="Yes",
                    Hamiltonian_SlaterKosterFiles_="Type2FileNames",
                    Hamiltonian_SlaterKosterFiles_Prefix=DFTB_PREFIX_3OB,
                    Hamiltonian_SlaterKosterFiles_Separator='"-"',
                    Hamiltonian_SlaterKosterFiles_Suffix='".skf"',
                    Hamiltonian_ThirdOrderFull="Yes",
                    Hamiltonian_HubbardDerivs_="",
                    Hamiltonian_HubbardDerivs_C="-0.1492",
                    Hamiltonian_HubbardDerivs_H="-0.1857",
                    **_dftb_common(charge))

    if method == "IPEA1":                             # xTB-IPEA via DFTB+
        from ase.calculators.dftb import Dftb
        return Dftb(label=os.path.join(workdir, "dftb"),
                    Hamiltonian_="xTB",
                    Hamiltonian_Method="IPEA1-xTB",
                    **_dftb_common(charge))

    if method == "g-xTB":
        Gxtb = load_gxtb_class()
        if Gxtb is None:
            raise SystemExit("ASE calculator for g-xTB not found; "
                             "set GXTB_ASE_MODULE to gxtb.py")
        return Gxtb(charge=charge, uhf=0,
                    command=GXTB_COMMAND if forces else GXTB_COMMAND_ENERGY,
                    label=os.path.join(workdir, "gxtb"))

    raise SystemExit(f"unknown method: {method}")


def available(method):
    """Check that the binary and the parameters needed are present."""
    if method in COMPOSITE:
        return os.path.isfile(ORCA_BIN) or shutil.which(ORCA_BIN) is not None
    if method in SEMI_EMPIRICAL:
        return shutil.which("mopac") is not None
    if method in ("DFTB2", "DFTB3", "IPEA1"):
        if shutil.which("dftb+") is None:
            return False
        if method == "DFTB2":
            return os.path.isdir(DFTB_PREFIX_MIO)
        if method == "DFTB3":
            return os.path.isdir(DFTB_PREFIX_3OB)
        return True
    if method == "g-xTB":
        return (shutil.which(GXTB_COMMAND.split()[0]) is not None
                and load_gxtb_class() is not None)
    return False


# --------------------------------------------------------------------- timing

def reset_calc(calc):
    """Forget the previous result, so that the next trial recomputes.

    Calculator.reset() only exists on the old family of ASE calculators (Dftb,
    MOPAC, Gxtb). ORCA goes through GenericFileIOCalculator, which carries the
    same two attributes without the method: calling calc.reset() there raises
    AttributeError, and t_sp as well as t_grad fall to zero for the six
    composite methods. So we do by hand what reset() does.
    """
    calc.atoms = None
    calc.results = {}


def timed(label, fn, repeat=1):
    """Time fn. With repeat > 1, returns the median of the trials.

    The primitives t_sp and t_grad last a few hundredths of a second for the
    fast methods, where the start-up of the program and the noise of the clock
    dominate. Repeating them and taking the median makes the measurement
    usable. The start-up is kept within the measurement, since every gradient
    of the Hessian pays it too.
    """
    times, status = [], "ok"
    for _ in range(repeat):
        t0 = time.perf_counter()
        try:
            fn()
        except Exception as exc:
            status = f"ECHEC({type(exc).__name__})"
            times.append(time.perf_counter() - t0)
            break
        times.append(time.perf_counter() - t0)
    ordered = "  trials: " + " ".join(f"{t:.3f}" for t in times)
    times.sort()
    dt = times[len(times) // 2]
    spread = (f"  [{times[0]:.3f}-{times[-1]:.3f}]" + ordered) if repeat > 1 else ""
    print(f"    {label:14} {dt:10.3f} s   {status}{spread}")
    return dt, status


def start_geometry(molecule, method):
    """Optimised neutral geometry of the method; failing that, another one.

    The cost of an evaluation does not depend on the geometry used, provided it
    is reasonable: a gradient at a geometry displaced by 0.01 A costs as much as
    at equilibrium. The choice of the structure is therefore of no consequence
    to the measurement.
    """
    here = os.path.dirname(os.path.abspath(__file__))
    dirs = ("structures/optimized", os.path.join(here, os.pardir, "geometries"))
    for d in dirs:
        g = f"{d}/{molecule}_{method}_neutral_bfgs.xyz"
        if os.path.isfile(g):
            return g
    import glob
    for d in dirs:
        c = sorted(glob.glob(f"{d}/{molecule}_*_neutral_bfgs.xyz"))
        if c:
            return c[0]
    print(f"    no geometry available for {molecule}")
    return None


def run(method, molecule, with_opt=False, repeat=5, no_grad=False, tag=""):
    print(f"--- {method}")
    if not available(method):
        print("    binaire ou parametres absents, methode ignoree")
        return None

    geom = start_geometry(molecule, method)
    if geom is None:
        return None

    # The working directory carries the --tag label: two concurrent streams
    # measuring the same method on the same molecule (the HF-3C controls of the
    # September campaign) would otherwise write the same orca.out and make each
    # other fail with an ASE ParseError.
    base = f"timing/{molecule}_{method}{'_' + tag if tag else ''}"
    natoms = len(read(geom))

    # 1. energy alone
    a = read(geom); a.center()
    a.calc = make_calc(method, base + "_sp")
    t_sp, s_sp = timed("energie", lambda: (reset_calc(a.calc),
                                           a.get_potential_energy()),
                       repeat=repeat)

    # 2. energy + gradient: the primitive repeated by optimisation and Hessian.
    #    Skipped with --no-grad; t_grad and the estimate of the Hessian are then
    #    nan and t_sp alone is measured.
    n_grad_hess = 6 * natoms + 1
    t_grad, s_grad, t_hess_est = float("nan"), "non mesure", float("nan")
    if not no_grad:
        b = read(geom); b.center()
        b.calc = make_calc(method, base + "_grad", forces=True)
        t_grad, s_grad = timed("energie+grad", lambda: (reset_calc(b.calc),
                                                        b.get_forces()),
                               repeat=repeat)

        # cost of the Hessian: estimated, never measured. ase.vibrations
        # proceeds by finite differences of the forces, that is 6 N + 1 gradient
        # evaluations; measuring it would only repeat t_grad, for several hours
        # per method.
        t_hess_est = n_grad_hess * t_grad
        print(f"    {'hessienne':14} {t_hess_est:10.1f} s   estime "
              f"({n_grad_hess} x t_grad)")

    t_opt, nsteps, s_opt = float("nan"), -1, "non mesure"
    if with_opt:
        d = read(geom); d.center()
        d.calc = make_calc(method, base + "_opt", forces=True)
        dyn = BFGS(d, logfile=base + "_opt.log")
        t_opt, s_opt = timed("optimisation", lambda: dyn.run(fmax=MAX_FORCE))
        nsteps = dyn.get_number_of_steps()

    return dict(method=method, natoms=natoms, t_sp=t_sp, t_grad=t_grad,
                t_hess_est=t_hess_est, t_opt=t_opt, nsteps=nsteps,
                status=";".join((s_sp, s_grad, s_opt)))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--repeat", type=int, default=5,
                   help="number of trials for t_sp and t_grad, median kept "
                        "(default 5; use 1 for the slow methods)")
    p.add_argument("--with-opt", action="store_true",
                   help="also measure a complete optimisation (indicative: the "
                        "number of steps depends on the path, not the method)")
    p.add_argument("--no-grad", action="store_true",
                   help="measure the energy alone: t_grad and the estimate of "
                        "the Hessian are nan")
    p.add_argument("--tag", default="",
                   help="label inserted in the name of the output file and of "
                        "the working directory, so that concurrent streams "
                        "cannot write the same files (the name is stamped to "
                        "the second)")
    p.add_argument("--molecule", default="coronene",
                   help="reference molecule (default: coronene)")
    p.add_argument("methods", nargs="*", default=None,
                   help="methods to measure (default: all)")
    args = p.parse_args()

    # one single thread for the linear-algebra libraries: the measurement is
    # made on one core, and a library taking several of its own accord would
    # distort the ratio between methods.
    for v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
        os.environ[v] = "1"

    methods = args.methods or ALL_METHODS
    out = f"timings_{args.molecule}_serial" \
          f"{'_' + args.tag if args.tag else ''}" \
          f"_{time.strftime('%Y%m%d-%H%M%S')}.txt"
    with open(out, "w") as f:
        f.write(f"# molecule={args.molecule} nprocs=1 "
                f"repeat={args.repeat} grad={'non' if args.no_grad else 'oui'} "
                f"— ordres de grandeur\n")
        f.write("\t".join(["Method", "Molecule", "NAtoms", "NProcs",
                           "t_energy_s", "t_energy_grad_s",
                           "t_hessian_estimated_s", "t_optimisation_s",
                           "n_steps", "Status"]) + "\n")

    print(f"Molecule : {args.molecule}   coeurs : 1   sortie : {out}\n")
    for m in methods:
        r = run(m, args.molecule, with_opt=args.with_opt, repeat=args.repeat,
                no_grad=args.no_grad, tag=args.tag)
        if r is None:
            continue
        with open(out, "a") as f:
            f.write(f"{r['method']}\t{args.molecule}\t{r['natoms']}\t1\t"
                    f"{r['t_sp']:.3f}\t{r['t_grad']:.3f}\t{r['t_hess_est']:.1f}\t"
                    f"{r['t_opt']:.2f}\t{r['nsteps']}\t"
                    f"{r['status']}\n")

    print(f"\nTermine. Resultats dans {out}")


if __name__ == "__main__":
    main()

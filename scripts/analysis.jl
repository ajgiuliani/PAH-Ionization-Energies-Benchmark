#!/usr/bin/env julia
# =====================================================================
# analysis.jl
#
# Regenerates every table and figure of
#
#   A. Giuliani and H. R. Hrodmarsson,
#   "Ionization Energies of Polycyclic Aromatic Hydrocarbons:
#    Benchmarking Semi Empirical, Tight Binding and Composite Methods",
#   ChemPhysChem.
#
# from the raw energies deposited in data/. Nothing outside this
# repository is required.
#
# Usage (from the repository root):
#
#     julia scripts/analysis.jl
#
# Output is written to output/ :
#
#     table1_IE1.tex         table2_IE2.tex        table3_IEd.tex
#     tableS2_ionization.tex tableS3_metrics.tex   tableS4_zpe.tex
#     tableS6_topology.tex
#     box.png                gToC.png
#
# Tables S7 and S8 (resolution of the metrics) are written by
# resampling.jl, and Figure S1 (cost against size) by
# figure_cost_vs_size.jl.
#
# The .tex files contain the rows of the present work only, and the
# numbers alone: the tabular preamble, the thin space of the thousands
# and the $-$ of the negative values belong to the manuscript. The
# literature rows of Tables 1 to 3 of the article are quoted from the
# cited references and are not reproduced here, no more than the
# xTB-GFN2 row of Table S4, which is read from the COMPAS database.
#
# Table S1 is not written: it is the content of data/results_*.txt,
# already deposited, and no operation is performed on it. Table S5
# likewise reproduces data/spin_contamination_HF-3C.txt as it stands,
# the two differences it reports having been taken by the script that
# produced that file.
#
# Dependencies: DataFrames, CSV, Statistics, StatsPlots, Plots
#     julia -e 'using Pkg; Pkg.add(["DataFrames","CSV","StatsPlots","Plots"])'
#
# Licence: MIT (see scripts/LICENSE). The data are CC BY 4.0.
# =====================================================================

using DataFrames, CSV, Statistics, StatsPlots, Plots
using Plots.PlotMeasures
using Printf

# Paths are resolved relative to this file, so the script may be run from
# anywhere. Everything it reads lives in data/ ; everything it writes goes
# to output/, which is created on first run.
const ROOT = normpath(joinpath(@__DIR__, ".."))
const OUT  = joinpath(ROOT, "output")
isdir(OUT) || mkpath(OUT)

# ---------------------------------------------------------------------
# 1. Benchmark set and experimental reference values
# ---------------------------------------------------------------------
# Experimental first ionization energies are the evaluated values of the
# NIST Chemistry WebBook. Second and double ionization energies come from
# the ion-spectroscopy and charge-stripping literature cited in the article.
# NaN marks the absence of an experimental determination.
#
# The second ionization energy of naphthalene (14.6 +/- 0.5 eV) is
# deliberately excluded: it is not a direct measurement but a value derived
# by Schroeder et al. from the data of Tobita et al., and its true
# uncertainty is much larger than quoted. See the article for the discussion.
#
# The double ionization energy of anthracene (20.1 +/- 0.2 eV, Hartman et al.,
# ref. in the article) IS available and is retained: the double-ionization
# benchmark therefore rests on six molecules, not five. (An earlier version of
# this array marked it NaN by mistake.)

const MOLECULES = ["naphthalene", "fluorene", "anthracene", "pyrene",
                   "tetracene", "coronene", "ovalene", "dicoronylene"]
const CARBON    = [10, 13, 14, 16, 18, 24, 32, 48]

const IE1_EXP = [8.144, 7.91, 7.439, 7.426, 6.97, 7.29, 6.71, NaN]
const IE2_EXP = [  NaN,  NaN,  11.7,  11.7,  10.9, 10.95, 10.0, 9.1]
const IED_EXP = [ 21.5, 21.0,  20.1,  19.8,  18.6,  18.5,  NaN, NaN]

# Quoted experimental uncertainties (1 sigma, in eV), aligned with the arrays
# above. NaN wherever the reference value itself is absent or excluded. These
# are the parenthesised uncertainties of the cited measurements; ovalene's first
# ionization energy carries no quoted uncertainty in the source, and 0.05 eV is
# assumed (it changes nothing at the level of two decimals). They are used only
# by the resampling analysis (scripts/resampling.jl).
const IE1_UNC = [0.001, 0.02, 0.006, 0.001, 0.05, 0.03, 0.05, NaN]
const IE2_UNC = [  NaN,  NaN,   0.1,   0.1,  0.1, 0.05,  0.1, 0.2]
const IED_UNC = [  0.2,  0.2,   0.2,   0.2,  0.2,  0.2,  NaN, NaN]

# Order used throughout: approximate methods first, then composite methods.
const METHODS = ["DFTB2", "DFTB3", "IPEA1", "PM6", "PM7", "RM1", "g-xTB",
                 "HF-3C", "B97-3C", "R2SCAN-3C", "PBEH-3C", "B3LYP-3C", "WB97X-3C"]

# Display names (the raw files use the internal labels above).
const PRETTY = Dict("IPEA1" => "xTB-IPEA", "WB97X-3C" => "wB97X-3C")
pretty(m) = get(PRETTY, m, m)

# Empirical shift prescribed by the authors of the xTB-IPEA Hamiltonian.
# It is the only empirical correction applied anywhere in this work.
const EMP_CORR = Dict(m => (m == "IPEA1" ? 4.8455 : 0.0) for m in METHODS)

# ---------------------------------------------------------------------
# 2. Raw energies
# ---------------------------------------------------------------------
function load_energies()
    d = Dict{Tuple{String,String},NamedTuple}()
    for f in ("results_semi-empirical.txt", "results_g-xTB.txt", "results_composite.txt")
        path = joinpath(ROOT, "data", f)
        isfile(path) || error("missing data file: $path")
        for r in CSV.File(path, delim = '\t')
            d[(String(r.Method), String(r.Molecule))] = (
                E_n_Rn = r.E_neutral_Rn,  E_c_Rc = r.E_cation_Rc,
                E_c_Rn = r.E_cation_Rn,   E_d_Rd = r.E_dication_Rd,
                E_d_Rc = r.E_dication_Rc, E_d_Rn = r.E_dication_Rn,
                Z_n    = r.EPZ_neutral,   Z_c    = r.EPZ_cation,
                Z_d    = r.EPZ_dication)
        end
    end
    return d
end

const E = load_energies()

# ---------------------------------------------------------------------
# 3. Ionization energies, DeltaSCF
# ---------------------------------------------------------------------
# First  : neutral  -> cation    at their own relaxed geometries
# Second : cation   -> dication
# Double : neutral  -> dication  (the empirical shift enters twice)
#
# Zero-point energies are added to each state when the ZPE-corrected value
# is requested. g-xTB carries no frequencies, so its ZPE entries are zero
# and its corrected and uncorrected values coincide, as stated in the article.

function ie(m, mol, which; zpe::Bool)
    e = E[(m, mol)]; c = EMP_CORR[m]
    z(a, b) = zpe ? (a + b) : a
    which === :first  && return z(e.E_c_Rc, e.Z_c) - z(e.E_n_Rn, e.Z_n) -   c
    which === :second && return z(e.E_d_Rd, e.Z_d) - z(e.E_c_Rc, e.Z_c) -   c
    which === :double && return z(e.E_d_Rd, e.Z_d) - z(e.E_n_Rn, e.Z_n) - 2c
    error("unknown ionization: $which")
end

# Vertical values (Supporting Information), all at the neutral or cation geometry.
function vie(m, mol, which)
    e = E[(m, mol)]; c = EMP_CORR[m]
    which === :first  && return e.E_c_Rn - e.E_n_Rn -  c
    which === :second && return e.E_d_Rc - e.E_c_Rc -  c
    which === :double && return e.E_d_Rn - e.E_n_Rn - 2c
    error("unknown ionization: $which")
end

# Molecules retained for each comparison: those with an experimental value.
subset(exp_vals) = [i for i in eachindex(MOLECULES) if !isnan(exp_vals[i])]

const IDX = Dict(:first  => subset(IE1_EXP),
                 :second => subset(IE2_EXP),
                 :double => subset(IED_EXP))
const EXP = Dict(:first => IE1_EXP, :second => IE2_EXP, :double => IED_EXP)
const UNC = Dict(:first => IE1_UNC, :second => IE2_UNC, :double => IED_UNC)

"Quoted experimental uncertainties (eV) of the retained molecules, in IDX order."
uncertainties(which) = [UNC[which][i] for i in IDX[which]]

"Signed errors of method `m` against experiment, ZPE-corrected."
errors(m, which) = [ie(m, MOLECULES[i], which; zpe = true) - EXP[which][i]
                    for i in IDX[which]]

# ---------------------------------------------------------------------
# 4. Statistical metrics (Table S3)
# ---------------------------------------------------------------------
bias(d)  = sum(d) / length(d)                                   # mean signed error
sigma(d) = sqrt(sum((d .- bias(d)).^2) / (length(d) - 1))       # scatter about the bias
mae(d)   = sum(abs.(d)) / length(d)
rmse(d)  = sqrt(sum(d .* d) / length(d))
maxad(d) = maximum(abs.(d))

function metrics_table()
    rows = DataFrame(method = String[], ionization = String[], MSE = Float64[],
                     sigma = Float64[], MAE = Float64[], RMSE = Float64[], MAX = Float64[])
    for (which, label) in ((:first, "IE1"), (:second, "IE2"), (:double, "IEd"))
        for m in METHODS
            d = errors(m, which)
            push!(rows, (pretty(m), label, bias(d), sigma(d), mae(d), rmse(d), maxad(d)))
        end
    end
    rows
end

# ---------------------------------------------------------------------
# 5. LaTeX output
# ---------------------------------------------------------------------
fmt(x) = isnan(x) ? "" : @sprintf("%.2f", x)

function write_ie_table(which, fname, caption)
    open(joinpath(OUT, fname), "w") do io
        println(io, "% ", caption)
        println(io, "% Present work only; literature rows are quoted in the article.")
        println(io, "% Columns: ", join(MOLECULES, " & "))
        for m in METHODS
            vals = [fmt(ie(m, mol, which; zpe = true)) for mol in MOLECULES]
            println(io, pretty(m), " & ", join(vals, " & "), " \\\\")
        end
    end
end

function write_metrics_table()
    df = metrics_table()
    open(joinpath(OUT, "tableS3_metrics.tex"), "w") do io
        println(io, "% Table S3 - statistical metrics, ZPE-corrected adiabatic values")
        println(io, "% Methods & MSE & MAE & RMSE & sigma")
        for label in ("IE1", "IE2", "IEd")
            println(io, "\\hline\n", label, " with ZPE & & & & \\\\\n\\hline")
            for r in eachrow(df[df.ionization .== label, :])
                @printf(io, "  %s & %.2f & %.2f & %.2f & %.2f \\\\\n",
                        r.method, r.MSE, r.MAE, r.RMSE, r.sigma)
            end
        end
    end
    df
end

# Table S2: the nine ionization energies of each molecule and method, one row
# per pair. It is the only table of the article to report the adiabatic values
# without the zero-point correction and the vertical values, and the only one to
# keep the internal method labels rather than the display names.
#
# Its numbers are printed to two decimals with the trailing zeros removed, as
# published: 8.00 is written 8 and 7.80 is written 7.8.
trim2(x) = String(rstrip(rstrip(@sprintf("%.2f", x), '0'), '.'))

function write_all_ie_table()
    open(joinpath(OUT, "tableS2_ionization.tex"), "w") do io
        println(io, "% Table S2 - adiabatic ionization energies with and without the zero-point")
        println(io, "% correction, and vertical ionization energies.")
        println(io, "% Molecule & Method & aIE1+ZPE & aIE1 & vIE1 & aIE2+ZPE & aIE2 & vIE2",
                    " & aIEd+ZPE & aIEd & vIEd")
        for mol in MOLECULES, m in METHODS
            vals = String[]
            for which in (:first, :second, :double)
                push!(vals, trim2(ie(m, mol, which; zpe = true)))
                push!(vals, trim2(ie(m, mol, which; zpe = false)))
                push!(vals, trim2(vie(m, mol, which)))
            end
            println(io, "  ", mol, " & ", m, " & ", join(vals, " & "), " \\\\")
        end
    end
end

# Table S4: zero-point energy of the neutral, and the correction DeltaZPE that
# enters the first adiabatic ionization energy. g-xTB is absent from the table
# for want of frequencies; its three ZPE columns are zero in the raw file and a
# row of zeros would be read as a result rather than as an absence.
#
# DeltaZPE is taken at full precision and rounded once, at printing. The
# published table does not do so throughout: its DFTB2, xTB-IPEA, PM6, PM7 and
# RM1 rows are the difference of the two ZPEs already rounded to three decimals,
# which shifts fifteen of the ninety-six values by 0.001 eV. The DFTB3 row and
# the six composite rows are at full precision, as here. The likely history is
# the correction of the DFTB3 zero-point energies, after which that row was
# recomputed and the four other approximate ones were not.
const ZPE_METHODS = filter(!=("g-xTB"), METHODS)

function write_zpe_table()
    open(joinpath(OUT, "tableS4_zpe.tex"), "w") do io
        println(io, "% Table S4 - zero-point energies of the neutral species and")
        println(io, "% DeltaZPE = ZPE(cation) - ZPE(neutral), in eV.")
        println(io, "% g-xTB is absent: the preliminary version used gives no frequencies.")
        println(io, "% The xTB-GFN2 row of the article comes from the COMPAS database and is")
        println(io, "% not derivable from data/ ; it is not written here.")
        println(io, "% Columns: ", join(MOLECULES, " & "))
        println(io, "% --- ZPE of the neutral species")
        for m in ZPE_METHODS
            vals = [@sprintf("%.3f", E[(m, mol)].Z_n) for mol in MOLECULES]
            println(io, pretty(m), " & ", join(vals, " & "), " \\\\")
        end
        println(io, "% --- DeltaZPE = ZPE(cation) - ZPE(neutral)")
        for m in ZPE_METHODS
            vals = [@sprintf("%.3f", E[(m, mol)].Z_c - E[(m, mol)].Z_n) for mol in MOLECULES]
            println(io, pretty(m), " & ", join(vals, " & "), " \\\\")
        end
    end
end

# Table S6: the two structural families, at comparable size ranges.
const ACENES        = ["naphthalene", "anthracene", "tetracene"]   # C10, C14, C18
const PERICONDENSED = ["pyrene", "coronene", "ovalene"]            # C16, C24, C32

function write_topology_table()
    dev(m, mol) = ie(m, mol, :first; zpe = true) - IE1_EXP[findfirst(==(mol), MOLECULES)]
    open(joinpath(OUT, "tableS6_topology.tex"), "w") do io
        println(io, "% Table S6 - deviation from experiment of the first adiabatic IE")
        println(io, "% linear acenes | pericondensed | drift over each subseries")
        for m in METHODS
            a = [dev(m, x) for x in ACENES]
            p = [dev(m, x) for x in PERICONDENSED]
            @printf(io, "%s & %.2f & %.2f & %.2f & & %.2f & %.2f & %.2f & & %.2f & %.2f \\\\\n",
                    pretty(m), a..., p..., a[end] - a[1], p[end] - p[1])
        end
    end
end

# ---------------------------------------------------------------------
# 6. Figure 2 - box plots of the errors
# ---------------------------------------------------------------------
# Whiskers span the full range (minimum to maximum): with five to seven
# molecules per box the Tukey 1.5 x IQR rule would flag ordinary points as
# outliers and hide them. Panels d) and e) carry a broken error axis so that
# the HF-3C deviations remain visible without compressing the other methods.

function figure_boxplots()
    d1 = Dict(m => errors(m, :first)  for m in METHODS)
    d2 = Dict(m => errors(m, :second) for m in METHODS)
    d3 = Dict(m => errors(m, :double) for m in METHODS)
    # Figure 2 keeps the internal method labels (IPEA1), unlike the heat map.
    labels = [m == "WB97X-3C" ? "ωB97X-3C" : m for m in METHODS]
    approx, comp = 1:7, 8:13

    function simple(dd, idx, ttl, loc; mirror = false)
        p = boxplot(xlabel = "error (eV)", permute = (:x, :y), title = ttl,
                    titlefont = font(10), titlelocation = loc, ymirror = mirror,
                    palette = mirror ? :darktest : :auto)
        for (i, m) in enumerate(METHODS[idx])
            boxplot!(p, [i], dd[m], label = false, outliers = false,
                     whisker_range = Inf, permute = (:x, :y))
        end
        vline!(p, [0], linestyle = :dash, c = :black, label = false)
        for (i, m) in enumerate(METHODS[idx])
            scatter!(p, [i], [mean(dd[m])], permute = (:x, :y), label = false, c = :black)
        end
        plot!(p, yticks = (1:length(idx), labels[idx]), ymirror = mirror)
        p
    end

    # broken-axis panel: `xl` is the window, `show_labels` the method names
    function broken(dd, xl, xtk; show_labels, xlab, ttl, yaxis = true)
        p = boxplot(xlabel = xlab, permute = (:x, :y), ymirror = true, xlims = xl,
                    title = ttl, titlefont = font(10), titlelocation = :right,
                    palette = :darktest, xticks = xtk, tickfont = font(7))
        for (i, m) in enumerate(METHODS[comp])
            boxplot!(p, [i], dd[m], label = false, outliers = false,
                     whisker_range = Inf, permute = (:x, :y))
        end
        for (i, m) in enumerate(METHODS[comp])
            scatter!(p, [i], [mean(dd[m])], permute = (:x, :y), label = false, c = :black)
        end
        if yaxis
            plot!(p, yticks = (1:length(comp),
                               show_labels ? labels[comp] : fill("", length(comp))),
                  ymirror = true)
        else
            plot!(p, yshowaxis = false, yticks = false, ygrid = false)
        end
        p
    end

    t(s) = plot(title = s, grid = false, showaxis = false, bottom_margin = -50Plots.px)

    p1 = simple(d1, approx, "a)", :left)
    p4a = broken(d1, (-2.25, -1.28), [-2.0, -1.5]; show_labels = false, xlab = "",
                 ttl = "d)", yaxis = false)
    p4b = broken(d1, (-0.72, 0.30), [-0.6, -0.3, 0.0]; show_labels = true,
                 xlab = "error (eV)", ttl = "")
    vline!(p4b, [0], linestyle = :dash, c = :black, label = false)
    annotate!(p4b, -0.72 + 0.010 * 1.02, 0.50, text("//", 11, :black, rotation = 72))

    p2 = simple(d2, approx, "b)", :left)
    p5a = broken(d2, (-0.55, 1.30), [-0.5, 0.0, 0.5, 1.0]; show_labels = false,
                 xlab = "error (eV)", ttl = "", yaxis = false)
    vline!(p5a, [0], linestyle = :dash, c = :black, label = false)
    p5b = broken(d2, (3.85, 4.20), [4.0]; show_labels = true, xlab = "", ttl = "e)")
    annotate!(p5b, 3.85 + 0.010 * 0.35, 0.50, text("//", 11, :black, rotation = 72))

    p3 = simple(d3, approx, "c)", :left)
    p6 = simple(d3, comp, "f)", :right; mirror = true)

    fig = plot(t("First IE"), p1, p4a, p4b, t("Second IE"), p2, p5a, p5b,
               t("Double IE"), p3, p6,
               layout = @layout([A{0.01h}; [B{0.50w} C{0.16w} D];
                                 A{0.01h}; [E{0.50w} F{0.38w} G];
                                 A{0.01h}; [H I]]),
               size = (600, 800), right_margin = 5mm, left_margin = 5mm)
    plot!(fig[3], right_margin = -2mm); plot!(fig[4], left_margin = -2mm)
    plot!(fig[7], right_margin = -2mm); plot!(fig[8], left_margin = -2mm)
    savefig(fig, joinpath(OUT, "box.png"))
end

# ---------------------------------------------------------------------
# 7. Graphical table of contents - MAE heat map
# ---------------------------------------------------------------------
# Sequential single-hue ramp, light to dark, spanning the full data range:
# magnitude is carried by lightness alone, so the map survives greyscale
# reproduction and colour-vision deficiency.

function figure_gtoc(df)
    order = reverse(METHODS)                       # row 1 of the heat map is the bottom
    lab   = [m == "WB97X-3C" ? "ωB97X-3C" : pretty(m) for m in order]
    mat = [only(df[(df.method .== pretty(order[i])) .&
                   (df.ionization .== ("IE1", "IE2", "IEd")[j]), :MAE])
           for i in eachindex(order), j in 1:3]
    n, vmax = length(order), maximum(mat)
    grad = cgrad([colorant"#cde2fb", colorant"#9ec5f4", colorant"#6da7ec",
                  colorant"#3987e5", colorant"#256abf", colorant"#184f95",
                  colorant"#0d366b"])

    # numeric axes: with string tick labels Plots offsets the annotations by
    # half a cell relative to the cells themselves
    p1 = heatmap(1:3, 1:n, mat, color = grad, clims = (0.0, vmax), colorbar = false,
                 aspect_ratio = 0.4, framestyle = :box, xmirror = true,
                 tickfontsize = 10, xticks = (1:3, ["IE₁", "IE₂", "IEd"]),
                 yticks = (1:n, lab))
    for i in 1:n, j in 1:3
        annotate!(p1, j, i, text(string(round(mat[i, j], digits = 2)), 8,
                                 mat[i, j] > 0.55 * vmax ? :white : :black, :center))
    end
    hline!(p1, [6.5], color = :gray40, linewidth = 1.5, linestyle = :dash, label = false)

    zbar = range(0, vmax, length = 200)
    p2 = heatmap([1], zbar, reshape(zbar, :, 1), color = grad, clims = (0.0, vmax),
                 colorbar = false, xticks = false, ymirror = true, ylabel = "MAE (eV)",
                 tickfontsize = 9, guidefontsize = 10, framestyle = :box,
                 yticks = ([0.0, 0.4, 0.8, 1.2, 1.6], ["0", "0.4", "0.8", "1.2", "1.6"]))
    savefig(plot(p1, p2, layout = grid(1, 2, widths = [0.88, 0.12]), size = (558, 560)),
            joinpath(OUT, "gToC.png"))
end

# ---------------------------------------------------------------------
# 8. Run
# ---------------------------------------------------------------------
function main()
    write_ie_table(:first,  "table1_IE1.tex", "Table 1 - first adiabatic IE, ZPE-corrected")
    write_ie_table(:second, "table2_IE2.tex", "Table 2 - second adiabatic IE, ZPE-corrected")
    write_ie_table(:double, "table3_IEd.tex", "Table 3 - double adiabatic IE, ZPE-corrected")
    write_all_ie_table()
    df = write_metrics_table()
    write_zpe_table()
    write_topology_table()
    figure_boxplots()
    figure_gtoc(df)

    println("\nWritten to ", OUT, " :")
    for f in sort(readdir(OUT)); println("  ", f); end
    println("\nMolecules per comparison: IE1 ", length(IDX[:first]),
            ", IE2 ", length(IDX[:second]), ", IEd ", length(IDX[:double]))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

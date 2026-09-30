#!/usr/bin/env julia
# =====================================================================
# figure_cost_vs_size.jl
#
# Cost of one single-point energy evaluation against the number of
# atoms, log-log, thirteen methods, four sizes:
#
#     coronene              C24H12    36 atoms
#     ovalene               C32H14    46 atoms
#     dicoronylene          C48H20    68 atoms
#     circumcircumcoronene  C96H24   120 atoms
#
# Reads data/timings_size_scaling.txt and writes
# output/figure_cost_vs_size.png (300 dpi) and .pdf.
#
# Usage (from anywhere):
#
#     julia scripts/figure_cost_vs_size.jl
#
# The times are cold single points (NoAutoStart) on one core; see the
# header of the data file. Colours are those of the other figures of
# the manuscript, one per method; circles for the fast methods, squares
# for the composite ones; in brackets after each label, the effective
# exponent alpha of the fit of ln t against ln N.
#
# Dependencies: Plots (GR backend)
#     julia -e 'using Pkg; Pkg.add("Plots")'
#
# Licence: MIT (see scripts/LICENSE). The data are CC BY 4.0.
# =====================================================================

using Plots, DelimitedFiles, Printf
gr()

# Paths are resolved relative to this file, as in analysis.jl, so the
# script may be run from anywhere.
const ROOT = normpath(joinpath(@__DIR__, ".."))
const OUT  = joinpath(ROOT, "output")
isdir(OUT) || mkpath(OUT)

const N = [36, 46, 68, 120]
const formulas = ["C₂₄H₁₂", "C₃₂H₁₄", "C₄₈H₂₀", "C₉₆H₂₄"]
# colours of the manuscript; the order is that of the labels
const color = Dict(
    "DFTB2" => "#009BFA", "DFTB3" => "#E36F47", "IPEA1" => "#3EA44D",
    "PM6" => "#C371D3", "PM7" => "#AC8E18", "RM1" => "#00AAAE", "g-xTB" => "#EE5E93",
    "HF-3C" => "#00008B", "B97-3C" => "#8A2BE2", "R2SCAN-3C" => "#008B8B",
    "PBEH-3C" => "#008000", "B3LYP-3C" => "#B3B300", "WB97X-3C" => "#FFA500")
const pretty = Dict("WB97X-3C" => "ωB97X-3C", "IPEA1" => "xTB-IPEA")
const fast = Set(["DFTB2", "DFTB3", "IPEA1", "PM6", "PM7", "RM1", "g-xTB"])

# ---- data --------------------------------------------------------------------
raw = readdlm(joinpath(ROOT, "data", "timings_size_scaling.txt"), '\t', String;
              comments = true, comment_char = '#')
hdr = raw[1, :]; rows = raw[2:end, :]
col(name) = findfirst(==(name), hdr)
methods = String.(rows[:, 1])
t = Dict(m => parse.(Float64, rows[i, 2:5]) for (i, m) in enumerate(methods))
alpha = Dict(m => parse(Float64, rows[i, col("alpha")]) for (i, m) in enumerate(methods))

# ---- figure ------------------------------------------------------------------
plt = plot(size=(640, 540), dpi=300, framestyle=:box, legend=false,
           xscale=:log10, yscale=:log10,
           xlims=(31, 345), ylims=(0.05, 1.2e5),
           xticks=(N, string.(N)), yticks=(10.0 .^ (-1:5), ["0.1", "1", "10", "100", "1000", "10⁴", "10⁵"]),
           xlabel="Number of atoms", ylabel="Single-point wall time (s)",
           minorgrid=false, gridalpha=0.25, gridstyle=:dot,
           guidefontsize=12, tickfontsize=10, left_margin=4Plots.mm, bottom_margin=3Plots.mm,
           right_margin=2Plots.mm, top_margin=2Plots.mm)

# formula above each size
for (n, f) in zip(N, formulas)
    annotate!(plt, n, 7.0e4, text(f, 9, :gray20, :center))
end

# curves: circle for the fast methods, square for the composite ones
order = sort(methods, by=m -> t[m][end])
for m in order
    plot!(plt, N, t[m], color=color[m], linewidth=2,
          marker=(m in fast ? :circle : :square), markersize=5,
          markerstrokecolor=:black, markerstrokewidth=0.6)
end

# labels on the right, spread out in log so that they do not overlap
mindist = 0.165                                   # in decades (height of a 9 pt label)
ys = [log10(t[m][end]) for m in order]
lab = copy(ys)
for _ in 1:200
    moved = false
    for i in 2:length(lab)
        if lab[i] - lab[i-1] < mindist
            d = (mindist - (lab[i] - lab[i-1])) / 2
            lab[i] += d; lab[i-1] -= d; moved = true
        end
    end
    moved || break
end
for (m, y0, y1) in zip(order, ys, lab)
    name = get(pretty, m, m)
    plot!(plt, [120 * 1.05, 120 * 1.18], [10^y0, 10^y1], color=:gray50, linewidth=0.5)
    annotate!(plt, 120 * 1.22, 10^y1, text(@sprintf("%s  (%.1f)", name, alpha[m]), 9, :black, :left))
end

savefig(plt, joinpath(OUT, "figure_cost_vs_size.png"))
savefig(plt, joinpath(OUT, "figure_cost_vs_size.pdf"))
println("written: ", joinpath(OUT, "figure_cost_vs_size.png"), " / .pdf")

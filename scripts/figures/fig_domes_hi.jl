# Figure res_domes_hi (thesis fig:domes_hi) — the generalized entropy at all four couplings in one
# uniform layout: one row per coupling, Re S2 with the c=1/2 curve of Eq. (G2) beside Im S2 with
# the pi*c/16 line. Each row spans that coupling's entropy window in Table 1.
# Reads: data/local/controls/seedens_p{0.0,0.1,0.3,0.5}_T*.jld2 (the single-vector seed ensembles)
# Run:   julia --project=. scripts/figures/fig_domes_hi.jl
# From a notebook: include this file, then call make_domes_hi() to get the plot object.

ROOT = normpath(joinpath(@__DIR__, "..", ".."))
isdefined(Main, :thesis_plot_theme!) || include(joinpath(ROOT, "src", "thesislib.jl"))
using JLD2, Printf, Statistics
thesis_plot_theme!()

wch(t, T) = (2T / pi) * sin(pi * t / T)
W(t, T) = log(wch(t, T))
# the c=1/2 profile of Eq. (G2): the offset and the amplitude of the w^(-1/2) correction are
# non-universal, so both are matched to each profile over its central 90%
function chord_pars(s2, ts, T)
    n = length(s2); lo = max(1, round(Int, n * 0.05) + 1); win = lo:(n + 1 - lo)
    w = wch.(ts, T)
    return [ones(length(win)) w[win] .^ (-0.5)] \ (real.(s2)[win] .- log.(w[win]) ./ 16)
end
Tlab(T) = T == round(T) ? string(Int(T)) : string(T)

# every rung of this coupling that has at least one physical seed
ent_arm(p) = Dict(T => r for T in 2.0:0.5:24.0
                  for r in (ensemble_profile(p, T),) if r !== nothing)
# times over which the profile still follows the chord line; past them the power method has
# left the physical branch and the profile is an artefact
const ENT_MAX = Dict(0.0 => 22.0, 0.1 => 13.0, 0.3 => 9.0, 0.5 => 5.0)   # windows of Table 1
# lower end: every rung entering a fit has v*T >~ 8, so the strip is wide enough in conformal units.
# Only p = 0 is cut by this, where v = 2 puts T = 2 and 3 at strips of 4 and 6 lattice units.
const ENT_MIN = Dict(0.0 => 4.0, 0.1 => 4.0, 0.3 => 2.0, 0.5 => 2.0)
const PLAB = Dict(0.0 => "0", 0.1 => "0.1", 0.3 => "0.3", 0.5 => "0.5")

# the times shown for one coupling: at most six, evenly spaced, so the legend stays clear
function dome_times(a, p)
    Ts = sort([T for T in keys(a) if ENT_MIN[p] <= T <= ENT_MAX[p]])
    length(Ts) > 6 && (Ts = Ts[round.(Int, range(1, length(Ts), length=6))])
    return Ts
end

# one row of the figure: Re S2 with the chord curves, Im S2 with the plateau line
# an empty panel tag drops the "(a)" prefix, for the standalone one-coupling figures
tag(t, plab) = isempty(t) ? "p = $plab" : "($t) p = $plab"

function dome_column(armd, Tshow, plab, panels)
    pal = cgrad(:viridis, length(Tshow), categorical=true)
    allre = vcat([real.(armd[T].s2) for T in Tshow]...)
    lo, hi = extrema(allre); pad = 0.08 * (hi - lo)
    pa = plot(xlabel="t/T", ylabel="Re S₂", legend=false,
              ylims=(lo - pad, hi + pad), title=tag(panels[1], plab))
    pb = plot(xlabel="t/T", ylabel="Im S₂", legend=:outerright, xticks=[0.0, 0.5, 1.0],
              title=tag(panels[2], plab))
    for (i, T) in enumerate(Tshow)
        s2 = armd[T].s2
        n = length(s2)
        ts = range(T / (n + 1), T - T / (n + 1), length=n)
        q = chord_pars(s2, ts, T)
        scatter!(pa, ts ./ T, real.(s2), color=pal[i], ms=3, msw=0, label="T=$(Tlab(T))")
        tt = range(0.08, 0.92, length=200)
        plot!(pa, tt, W.(tt .* T, T) ./ 16 .+ q[1] .+ q[2] .* wch.(tt .* T, T) .^ (-0.5),
              color=:black, lw=1.0, ls=:dash, alpha=0.85, label="")
        scatter!(pb, ts ./ T, imag.(s2), color=pal[i], ms=3, msw=0, label="T=$(Tlab(T))")
    end
    hline!(pb, [pi / 32], color=:black, lw=2.0, ls=:dash, label="πc/16")
    return pa, pb
end

function make_domes_hi()
    hi = [(0.0, ("a", "b")), (0.1, ("c", "d")), (0.3, ("e", "f")), (0.5, ("g", "h"))]
    rows = Plots.Plot[]
    for (p, tags) in hi
        a = ent_arm(p)
        pa, pb = dome_column(a, dome_times(a, p), PLAB[p], tags)
        push!(rows, pa); push!(rows, pb)
    end

    return plot(rows..., layout=@layout([a b{0.545w}; c d{0.545w}; e f{0.545w}; g h{0.545w}]),
                size=thesis_size(1.0; aspect=1.25), margin=3Plots.mm,
                bottom_margin=6Plots.mm, left_margin=7Plots.mm)
end

# One coupling on its own, in the same format as a row of the full figure. The defence deck shows
# the four couplings as separate images rather than as one grid, and building them here keeps them
# in step with Figure fig:domes_hi instead of relying on crops of it.
function make_dome_single(p)
    a = ent_arm(p)
    pa, pb = dome_column(a, dome_times(a, p), PLAB[p], ("", ""))
    return plot(pa, pb, layout=@layout([a b{0.545w}]),
                size=thesis_size(0.95; aspect=0.40), margin=2Plots.mm,
                bottom_margin=4Plots.mm, left_margin=5Plots.mm)
end

if abspath(PROGRAM_FILE) == @__FILE__
    save_thesis_figure(make_domes_hi(), "res_domes_hi")
    for (p, tag) in ((0.0, "p0"), (0.1, "p01"), (0.3, "p03"), (0.5, "p05"))
        save_thesis_figure(make_dome_single(p), "res_domes_$tag")
    end
end

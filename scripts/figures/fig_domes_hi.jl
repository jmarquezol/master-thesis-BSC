# Figure res_domes_hi (thesis fig:domes_hi) — the generalized entropy at the couplings the main
# text does not show (p = 0, 0.3, 0.5), in the layout of fig:domes: one row per coupling, Re S2
# with the c=1/2 chord curves beside Im S2 with the pi*c/16 line.
# Reads: data/local/controls/seedens_p{0.0,0.3,0.5}_T*.jld2 (the single-vector seed ensembles)
# Run:   julia --project=. scripts/figures/fig_domes_hi.jl
# From a notebook: include this file, then call make_domes_hi() to get the plot object.

ROOT = normpath(joinpath(@__DIR__, "..", ".."))
isdefined(Main, :thesis_plot_theme!) || include(joinpath(ROOT, "src", "thesislib.jl"))
using JLD2, Printf, Statistics
thesis_plot_theme!()

function make_domes_hi()
    W(t, T) = log((2T / pi) * sin(pi * t / T))
    Tlab(T) = T == round(T) ? string(Int(T)) : string(T)

    # every rung of this coupling that has at least one physical seed
    ent_arm(p) = Dict(T => r for T in 2.0:0.5:24.0
                      for r in (ensemble_profile(p, T),) if r !== nothing)
    # times over which the profile still follows the chord line; past them the power method has
    # left the physical branch and the profile is an artefact
    ent_max = Dict(0.0 => 24.0, 0.3 => 8.0, 0.5 => 5.0)   # the entropy windows of Table 1

    # one row of the figure: Re S2 with the chord curves, Im S2 with the plateau line
    function dome_column(armd, Tshow, plab, panels)
        pal = cgrad(:viridis, length(Tshow), categorical=true)
        pa = plot(xlabel="t/T", ylabel="Re S₂", legend=false,
                  title="($(panels[1])) p = $plab")
        pb = plot(xlabel="t/T", ylabel="Im S₂", legend=:outerright, xticks=[0.0, 0.5, 1.0],
                  title="($(panels[2])) p = $plab")
        for (i, T) in enumerate(Tshow)
            s2 = armd[T].s2
            n = length(s2)
            ts = range(T / (n + 1), T - T / (n + 1), length=n)
            bulk = (n ÷ 4):(3n ÷ 4)
            # s0 is non-universal, so match it to the data and compare shapes only
            s0 = mean(real.(s2)[bulk] .- (0.5 / 8) .* W.(ts[bulk], T))
            scatter!(pa, ts ./ T, real.(s2), color=pal[i], ms=3, msw=0, label="T=$(Tlab(T))")
            tt = range(0.04, 0.96, length=200)
            plot!(pa, tt, (0.5 / 8) .* W.(tt .* T, T) .+ s0, color=pal[i], lw=1.2, ls=:dash, label="")
            scatter!(pb, ts ./ T, imag.(s2), color=pal[i], ms=3, msw=0, label="T=$(Tlab(T))")
        end
        hline!(pb, [pi / 32], color=:black, lw=2.0, ls=:dash, label="πc/16")
        return pa, pb
    end

    hi = [(0.0, "0", ("a", "b")), (0.3, "0.3", ("c", "d")), (0.5, "0.5", ("e", "f"))]
    rows = Plots.Plot[]
    for (p, plab, tags) in hi
        a = ent_arm(p)
        Ts = sort([T for T in keys(a) if T <= ent_max[p]])
        # at most six times per row, evenly spaced, so the legend stays clear of the axis
        length(Ts) > 6 && (Ts = Ts[round.(Int, range(1, length(Ts), length=6))])
        pa, pb = dome_column(a, Ts, plab, tags)
        push!(rows, pa)
        push!(rows, pb)
    end

    return plot(rows..., layout=@layout([a b{0.545w}; c d{0.545w}; e f{0.545w}]),
                size=thesis_size(1.0; aspect=0.95), margin=3Plots.mm,
                bottom_margin=6Plots.mm, left_margin=7Plots.mm)
end

if abspath(PROGRAM_FILE) == @__FILE__
    save_thesis_figure(make_domes_hi(), "res_domes_hi")
end

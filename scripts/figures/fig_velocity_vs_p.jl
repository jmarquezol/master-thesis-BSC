# Figure velocity_vs_p (thesis fig:velocity) — the sound velocity against the coupling: the
# ring values at every size N=10..18, with the N -> infinity extrapolation in black diamonds.
# These extrapolated values are the velocities used throughout the thesis.
# Reads: data/local/nb4_velocity_sizes.jld2, data/local/alcaraz_velocity.jld2
# Run:   julia --project=. scripts/figures/fig_velocity_vs_p.jl
# From a notebook: include this file, then call make_velocity_vs_p() to get the plot object.

ROOT = normpath(joinpath(@__DIR__, "..", ".."))
isdefined(Main, :thesis_plot_theme!) || include(joinpath(ROOT, "src", "thesislib.jl"))
using JLD2, Printf
thesis_plot_theme!()

function make_velocity_vs_p()
    N_ladder = (10, 12, 14, 16, 18)
    p_production = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0, 1.5]

    vN = load(joinpath(ROOT, "data", "local", "nb4_velocity_sizes.jld2"), "v")
    v_infinity = load(joinpath(ROOT, "data", "local", "alcaraz_velocity.jld2"), "v")

    figv = plot(size=thesis_size(0.92; aspect=0.45), xlabel="p", ylabel="v(p)",
                legend=:outerright, framestyle=:box)
    palette_sizes = cgrad(:viridis, length(N_ladder), categorical=true)

    for (i, N) in enumerate(N_ladder)
        scatter!(figv, p_production, [vN[(p, N)] for p in p_production],
                 color=palette_sizes[i], ms=4, msw=0, label="N = $N")
    end

    scatter!(figv, p_production, [v_infinity[p] for p in p_production],
             color=:black, marker=:diamond, ms=6, msw=0, label="extrapolated")
    hline!(figv, [2.0], color=:gray, ls=:dot, lw=1.2, label="exact, v = 2")

    @printf("v_inf(p): %s\n",
            join([@sprintf("v(%.1f)=%.3f", p, v_infinity[p]) for p in p_production], "  "))
    return figv
end

if abspath(PROGRAM_FILE) == @__FILE__
    save_thesis_figure(make_velocity_vs_p(), "velocity_vs_p")
end

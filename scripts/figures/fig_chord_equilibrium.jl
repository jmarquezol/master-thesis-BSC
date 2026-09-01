# Figure cft_chord_equilibrium (thesis fig:cft_chord) — equilibrium ground-state entanglement
# profiles at N=400 against the chord variable W(l,L), one curve per coupling, each with its
# fitted central charge in the legend. The profiles collapse on the Ising line c=1/2.
# Reads: data/local/nb4_chord_N400.jld2  (heavy DMRG, produced by notebooks/01)
# Run:   julia --project=. scripts/figures/fig_chord_equilibrium.jl
# From a notebook: include this file, then call make_chord_equilibrium() to get the plot object.

ROOT = normpath(joinpath(@__DIR__, "..", ".."))
isdefined(Main, :thesis_plot_theme!) || include(joinpath(ROOT, "src", "thesislib.jl"))
using JLD2, Printf, LsqFit
thesis_plot_theme!()

function make_chord_equilibrium()
    N = 400
    bonds = collect(1:N-1)
    p_chord = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0, 1.5]   # the full range this work studies
    chordmodel(x, q) = @. (q[1] / 6) * x + q[2]           # q[1] = c, q[2] = g

    cache = joinpath(ROOT, "data", "local", "nb4_chord_N$(N).jld2")
    isfile(cache) || error("missing $cache — the N=400 DMRG profiles; see notebooks/01 to regenerate")
    chord_data = load(cache, "data")

    figX = plot(xlabel="W(l, L)", ylabel="S_vN(x)", legend=:outerright,
                framestyle=:box, size=thesis_size(1.0; aspect=0.413))
    palette_chord = cgrad(:viridis, length(p_chord), categorical=true)

    for (i, p) in enumerate(p_chord)
        S    = chord_data[p]
        xv   = log.((2N / pi) .* sin.(pi .* bonds ./ N))
        bulk = findall(b -> N / 4 <= b <= 3N / 4, bonds)     # fit the bulk, skip the boundaries

        fit = curve_fit(chordmodel, xv[bulk], S[bulk], [0.5, 1.0])
        c   = fit.param[1]

        @printf("p=%.2f : c = %.4f\n", p, c)
        # one legend entry per coupling, carrying its fitted central charge
        scatter!(figX, xv, S; label="p=$p ($(round(c, digits=3)))", color=palette_chord[i], ms=3, alpha=0.55)
        plot!(figX, xv[bulk], chordmodel(xv[bulk], fit.param); label="", color=palette_chord[i], lw=2)
    end

    # spread of the profiles across p at the first and the central cut, both quoted in the thesis
    first_cut   = [chord_data[p][1] for p in p_chord]
    central_cut = [chord_data[p][div(N, 2)] for p in p_chord]
    @printf("spread across p: %.4f at the first cut, %.4f at the central cut\n",
            maximum(first_cut) - minimum(first_cut), maximum(central_cut) - minimum(central_cut))

    return figX
end

if abspath(PROGRAM_FILE) == @__FILE__
    save_thesis_figure(make_chord_equilibrium(), "cft_chord_equilibrium")
end

# Figure res_domes (thesis fig:domes) — the generalized entropy at p=0.1, real part beside
# imaginary part. Profiles come from the corrected-column cluster arm; the c=1/2 chord curve is
# drawn with the non-universal constant matched to each profile, and the Im panel is unzoomed
# with pi*c/16 marked.
# Reads: data/local/controls/seedens_p0.1_T*.jld2 (the single-vector seed ensembles)
# Run:   julia --project=. scripts/figures/fig_domes.jl
# From a notebook: include this file, then call make_domes() to get the plot object.

ROOT = normpath(joinpath(@__DIR__, "..", ".."))
isdefined(Main, :thesis_plot_theme!) || include(joinpath(ROOT, "src", "thesislib.jl"))
using JLD2, Printf, Statistics
thesis_plot_theme!()

function make_domes()
    W(t, T) = log((2T / pi) * sin(pi * t / T))
    plateau(s2) = mean(imag.(s2)[max(1, length(s2) ÷ 2 - 1):(length(s2) ÷ 2 + 2)])

    # Times that resolve: the plateau falls smoothly and the seed spread stays small up to T=13.
    # Beyond it the ensembles thin out and the plateau stops decreasing, so those rungs are left out.
    Tshow = [6.0, 8.0, 10.0, 12.0, 13.0]
    pal = cgrad(:viridis, length(Tshow), categorical=true)

    pa = plot(xlabel="t/T", ylabel="Re S₂", legend=false, title="(a)")
    pb = plot(xlabel="t/T", ylabel="Im S₂", legend=:outerright, title="(b)")
    for (i, T) in enumerate(Tshow)
        s2 = ensemble_profile(0.1, T).s2
        n = length(s2)
        ts = range(T / (n + 1), T - T / (n + 1), length=n)
        bulk = (n ÷ 4):(3n ÷ 4)
        s0 = mean(real.(s2)[bulk] .- (0.5 / 8) .* W.(ts[bulk], T))
        scatter!(pa, ts ./ T, real.(s2), color=pal[i], ms=3, msw=0, label="")
        tt = range(0.04, 0.96, length=200)
        plot!(pa, tt, (0.5 / 8) .* W.(tt .* T, T) .+ s0, color=pal[i], lw=1.2, ls=:dash, label="")
        scatter!(pb, ts ./ T, imag.(s2), color=pal[i], ms=3, msw=0, label="T=$(Int(T))")
        @printf("T=%2d  plateau = %.4f\n", Int(T), plateau(s2))
    end
    hline!(pb, [pi / 32], color=:black, lw=1.2, ls=:dot, label="πc/16")
    @printf("target πc/16 at c=1/2 = %.4f\n", pi / 32)

    return plot(pa, pb, layout=@layout([a b{0.589w}]),
                size=thesis_size(0.95; aspect=0.36), margin=2Plots.mm)
end

if abspath(PROGRAM_FILE) == @__FILE__
    save_thesis_figure(make_domes(), "res_domes")
end

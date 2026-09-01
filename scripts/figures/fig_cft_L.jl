# Figure cft_L (thesis fig:cft_L) — equilibrium finite-size scaling of the ground-state
# entanglement entropy with the corrected fits: (a) S(L/2) against ln L, (b) the subleading
# correction isolated. Hue is the estimator (S1 blue, S2 red), shade is the coupling.
# Reads: data/local/nb4_fss.jld2
# Run:   julia --project=. scripts/figures/fig_cft_L.jl
# From a notebook: include this file, then call make_cft_L() to get the plot object.

ROOT = normpath(joinpath(@__DIR__, "..", ".."))
isdefined(Main, :thesis_plot_theme!) || include(joinpath(ROOT, "src", "thesislib.jl"))
using JLD2, Printf
thesis_plot_theme!()

function make_cft_L()
    N_values = [40, 60, 80, 100, 120, 140, 160, 180, 200]
    fss_data = load(joinpath(ROOT, "data", "local", "nb4_fss.jld2"), "data")

    # S = m ln L + k + A L^(-x/n), linear in (m,k,A) at fixed x; x = 1 throughout
    function fit_corrected(S, n; x=1.0)
        L = Float64.(N_values)
        X = hcat(log.(L), ones(length(L)), L .^ (-x / n))
        return X \ S
    end

    figA = plot(xlabel="ln(L)", ylabel="S(L/2)", legend=false, xticks=[3.8, 4.2, 4.6, 5.0],
                framestyle=:box, title="(a)", titlelocation=:left, titlefontsize=11)
    figB = plot(xlabel="L", ylabel="S - (m ln L + k)", legend=false, xticks=[50, 100, 150, 200],
                framestyle=:box, title="(b)", titlelocation=:left, titlefontsize=11)
    # third panel carries the legend alone, so (a) and (b) keep identical plotting areas
    figLeg = plot(framestyle=:none, showaxis=false, grid=false, ticks=nothing, legend=:left,
                  foreground_color_legend=nothing, background_color_legend=nothing,
                  left_margin=6Plots.mm)

    # hue is the estimator, shade is the coupling: S1 in two blues, S2 in two reds
    styles = [(0.1, 1, 6.0, :deepskyblue, :circle, "S₁"),
              (0.5, 1, 6.0, :navy,        :square, "S₁"),
              (0.1, 2, 8.0, :salmon,      :circle, "S₂"),
              (0.5, 2, 8.0, :darkred,     :square, "S₂")]

    Lg = range(minimum(N_values), maximum(N_values); length=200)
    for (p, n, pref, col, mk, sym) in styles
        S = n == 1 ? collect(fss_data[p].vn) : collect(fss_data[p].r2)
        q = fit_corrected(S, n)
        c = pref * q[1]
        @printf("p=%.1f  n=%d :  c = %.4f   A = %+.3f\n", p, n, c, q[3])

        scatter!(figA, log.(N_values), S; label="", color=col, marker=mk, ms=4.5, msw=0)
        scatter!(figLeg, [NaN], [NaN]; label="  $sym, p=$p", color=col, marker=mk, ms=4.5, msw=0)
        plot!(figA, log.(Lg), q[1] .* log.(Lg) .+ q[2] .+ q[3] .* Lg .^ (-1 / n);
              label="", color=col, lw=2)

        scatter!(figB, N_values, S .- (q[1] .* log.(N_values) .+ q[2]);
                 label="", color=col, marker=mk, ms=4.5, msw=0)
        plot!(figB, Lg, q[3] .* Lg .^ (-1 / n); label="", color=col, lw=2)
    end

    return plot(figA, figB, figLeg; layout=@layout([a b c{0.20w}]),
                size=thesis_size(1.0; aspect=0.407))
end

if abspath(PROGRAM_FILE) == @__FILE__
    save_thesis_figure(make_cft_L(), "cft_L")
end

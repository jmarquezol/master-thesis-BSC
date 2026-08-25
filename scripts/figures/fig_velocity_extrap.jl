# Figure velocity_extrapolation (thesis fig:velocity_extrap) — two panels:
#   (a) every coupling, v(N) against 1/N^2, extrapolated leftwards to the intercept at N -> inf;
#   (b) the p = 0 control alone, where the intercept must land on the exact v = 2.
# Reads: data/local/nb4_velocity_sizes.jld2
# Run:   julia --project=. scripts/figures/fig_velocity_extrap.jl
# From a notebook: include this file, then call make_velocity_extrap() to get the plot object.

ROOT = normpath(joinpath(@__DIR__, "..", ".."))
isdefined(Main, :thesis_plot_theme!) || include(joinpath(ROOT, "src", "thesislib.jl"))
using JLD2, Printf
thesis_plot_theme!()

function make_velocity_extrap()
    N_ladder = (10, 12, 14, 16, 18)
    p_production = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0, 1.5]
    vN = load(joinpath(ROOT, "data", "local", "nb4_velocity_sizes.jld2"), "v")

    # v(N) = v_inf + C/N^2, linear in both parameters
    function fit_inv_sq(p)
        x = [1.0 / N^2 for N in N_ladder]
        y = [vN[(p, N)] for N in N_ladder]
        q = hcat(ones(length(x)), x) \ y
        return q[1], q[2]                      # v_inf, C
    end

    pal = cgrad(:viridis, length(p_production), categorical=true)
    xline = range(0, 1.05e-2, length=50)

    figA = plot(xlabel="10³ / N²", ylabel="v", legend=false, framestyle=:box,
                title="(a)", titlelocation=:left, titlefontsize=11)
    figB = plot(xlabel="10³ / N²", ylabel="v", legend=false, framestyle=:box,
                title="(b)  p = 0 control", titlelocation=:left, titlefontsize=11)
    figLeg = plot(framestyle=:none, showaxis=false, grid=false, ticks=nothing, legend=:left,
                  foreground_color_legend=nothing, background_color_legend=nothing,
                  left_margin=6Plots.mm)

    for (i, p) in enumerate(p_production)
        vinf, C = fit_inv_sq(p)
        x = [1000.0 / N^2 for N in N_ladder]
        y = [vN[(p, N)] for N in N_ladder]
        scatter!(figA, x, y; color=pal[i], ms=4, msw=0, label="")
        plot!(figA, 1000 .* xline, vinf .+ C .* xline; color=pal[i], lw=1.5, label="")
        scatter!(figA, [0.0], [vinf]; color=pal[i], marker=:diamond, ms=6, msw=0, label="")
        scatter!(figLeg, [NaN], [NaN]; label="  p = $p", color=pal[i], ms=4, msw=0)

        if p == 0.0                                  # the control, on its own scale
            scatter!(figB, x, y; color=pal[i], ms=5, msw=0, label="")
            plot!(figB, 1000 .* xline, vinf .+ C .* xline; color=pal[i], lw=1.8, label="")
            scatter!(figB, [0.0], [vinf]; color=pal[i], marker=:diamond, ms=8, msw=0, label="")
            hline!(figB, [2.0]; color=:black, ls=:dash, lw=1.4, label="")
            @printf("p=0: intercept %.5f, exact 2, error %.1e\n", vinf, abs(2 - vinf))
        end
    end

    # the intercept is the thermodynamic limit: label it where there is room, in the control panel
    annotate!(figB, 0.35, 1.9955, text("N → ∞", 9, :left))

    return plot(figA, figB, figLeg; layout=@layout([a b c{0.16w}]),
                size=thesis_size(1.0; aspect=0.40))
end

if abspath(PROGRAM_FILE) == @__FILE__
    save_thesis_figure(make_velocity_extrap(), "velocity_extrapolation")
end

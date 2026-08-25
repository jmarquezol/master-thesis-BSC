# Figure res_beta0 (thesis fig:beta0) — the fitted phase coefficient against the cooling
# regulator: it falls linearly with beta0 and extrapolates to the CFT prediction as beta0 -> 0.
# Reads: data/cluster/sweep_beta_p0.0.jld2
# Run:   julia --project=. scripts/figures/fig_beta0.jl
# From a notebook: include this file, then call make_beta0() to get the plot object.

ROOT = normpath(joinpath(@__DIR__, "..", ".."))
isdefined(Main, :thesis_plot_theme!) || include(joinpath(ROOT, "src", "thesislib.jl"))
using JLD2, Printf, Statistics, LsqFit
thesis_plot_theme!()

function make_beta0()
    ph(z) = angle(-z)
    dphw(a, b) = mod(a - b + pi, 2pi) - pi
    function unwrapped(mus)
        p = [ph(m) for m in mus]
        for i in 2:length(p); p[i] = p[i-1] + dphw(p[i], p[i-1]); end
        return p
    end
    eq3(x, q) = @. q[1] - pi / x + q[2] / x^2
    CPRED = -pi / 2 * (0.5 / 24)               # CFT prediction for the 1/T^2 coefficient at c=1/2

    b = load(joinpath(ROOT, "data", "cluster", "sweep_beta_p0.0.jld2"), "done")
    nb_of(k) = parse(Int, match(r"nb(\d+)", k[1]).captures[1])
    b0s = Float64[]; ratios = Float64[]
    for nb in sort(unique(nb_of.(keys(b))))
        ks = sort([k for k in keys(b) if nb_of(k) == nb], by=k -> k[2])
        Tb = [k[2] for k in ks]; length(Tb) < 5 && continue
        y = unwrapped([b[k].theta_phys for k in ks]) ./ Tb
        f = curve_fit(eq3, Tb, y, [0.1, -0.01])
        push!(b0s, nb * 0.05); push!(ratios, f.param[2] / CPRED)
    end
    lin(x, q) = @. q[1] + q[2] * x
    fit = curve_fit(lin, b0s, ratios, [1.0, -1.0])
    @printf("intercept at beta0=0: %.3f   slope %.3f\n", fit.param[1], fit.param[2])

    xs = range(0, maximum(b0s) * 1.05, length=50)
    pl = plot(xlabel="\$\\beta_0\$", ylabel="\$C_{\\mathrm{fit}}\\,/\\,C_{\\mathrm{CFT}}\$",
              legend=:topright, legendfontsize=7, size=thesis_size(0.62; aspect=0.72),
              margin=2Plots.mm, ylims=(0.4, 1.12), xlims=(-0.02, 0.86))
    hline!(pl, [1.0], color=:black, ls=:dash, lw=0.9, label="")
    plot!(pl, xs, lin(collect(xs), fit.param), color=:grey45, lw=1.4, label="linear fit")
    scatter!(pl, b0s, ratios, color=:dodgerblue, ms=4, msw=0, label="measured")
    scatter!(pl, [0.0], [fit.param[1]], color=:crimson, ms=5, msw=0, marker=:diamond,
             label="extrapolated, \$\\beta_0\\to0\$")
    scatter!(pl, [0.2], [ratios[findfirst(==(0.2), b0s)]], color=:seagreen, ms=6, msw=0,
             marker=:star5, label="production, \$\\beta_0=0.2\$")
    return pl
end

if abspath(PROGRAM_FILE) == @__FILE__
    save_thesis_figure(make_beta0(), "res_beta0")
end

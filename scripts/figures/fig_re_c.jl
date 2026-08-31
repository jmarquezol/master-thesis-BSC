# Figure res_re_c (thesis fig:re_c) — the central charge from the real part of the temporal
# entropy: one profile with the two fits it admits. The corrected c=1/2 form of Eq. (G2) describes
# the data, while a plain logarithm with c free follows the middle of the strip and leaves it in
# the wings, returning a c far above 1/2.
# Reads: data/local/controls/seedens_p0.0_T*.jld2 (the single-vector seed ensembles)
# Run:   julia --project=. scripts/figures/fig_re_c.jl
# From a notebook: include this file, then call make_re_c() to get the plot object.

ROOT = normpath(joinpath(@__DIR__, "..", ".."))
isdefined(Main, :thesis_plot_theme!) || include(joinpath(ROOT, "src", "thesislib.jl"))
using JLD2, Printf, Statistics
thesis_plot_theme!()

function make_re_c()
    wch(t, T) = (2T / pi) * sin(pi * t / T)
    cuts(n, T) = (dt = T / (n + 1); collect(dt:dt:(n * dt)))
    window(n, frac) = (lo = max(1, round(Int, n * (1 - frac) / 2) + 1); lo:(n + 1 - lo))
    lsq(X, y) = X \ y

    d = Dict(T => r for T in 2.0:0.5:24.0
             for r in (ensemble_profile(0.0, T),) if r !== nothing)
    Ts = sort(collect(keys(d)))

    function fits(s2, T)
        n = length(s2); ts = cuts(n, T)
        W = log.(wch.(ts, T)); cv = wch.(ts, T) .^ (-0.5); re = real.(s2)
        m50 = window(n, 0.5); c90 = window(n, 0.9)
        q1 = lsq([ones(length(m50)) W[m50]], re[m50])
        q3 = lsq([ones(length(c90)) W[c90] cv[c90]], re[c90])
        q2 = lsq([ones(length(c90)) cv[c90]], re[c90] .- W[c90] ./ 16)
        return (ts=ts, W=W, cv=cv, re=re, m50=m50, c90=c90,
                c_unc=8q1[2], s_unc=q1[1], c_corr=8q3[2], q3=q3, s_pin=q2[1], a_pin=q2[2])
    end

    Tshow = 14.0
    f = fits(d[Tshow].s2, Tshow)
    xs = f.ts ./ Tshow
    pa = plot(xlabel="\$t/T\$", ylabel="\$\\mathrm{Re}\\,S_2\$", title="\$T=$(Int(Tshow))\$", titlefontsize=11, legend=:bottom, legendfontsize=8,
              ylims=(0.05, 0.42))
    plot!(pa, xs, f.re, seriestype=:scatter, ms=2.5, msw=0, color=:grey40, label="data")
    plot!(pa, xs[f.c90], f.s_pin .+ f.W[f.c90] ./ 16 .+ f.a_pin .* f.cv[f.c90], color=:dodgerblue, lw=1.8,
          label="\$c=1/2\$ with correction")
    # fitted on the middle half, drawn across the wider window: it tracks the data where it was
    # fitted and leaves it in the wings
    plot!(pa, xs[f.c90], f.s_unc .+ (f.c_unc / 8) .* f.W[f.c90], color=:crimson, lw=1.6, ls=:dash,
          label="logarithm only, \$c=$(round(f.c_unc, digits=2))\$")

    cu = [fits(d[T].s2, T).c_unc for T in Ts]
    cc = [fits(d[T].s2, T).c_corr for T in Ts]
    @printf("T=%g: c_unc=%.3f  c_corr=%.3f  a_pin=%.4f\n", Tshow, f.c_unc, f.c_corr, f.a_pin)
    @printf("c_unc range %.3f..%.3f ; c_corr range %.3f..%.3f\n",
            minimum(cu), maximum(cu), minimum(cc), maximum(cc))

    # thesis_plot_theme! sets bottom_margin = 10mm globally, which a generic `margin` does not
    # override; unset here it leaves a blank band between the axis label and the caption
    plot!(pa, size=thesis_size(0.55; aspect=0.62), margin=2Plots.mm, bottom_margin=2Plots.mm)
    return pa
end

if abspath(PROGRAM_FILE) == @__FILE__
    save_thesis_figure(make_re_c(), "res_re_c")
end

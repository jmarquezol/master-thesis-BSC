# Figure res_conv (thesis fig:conv) — how the single-vector iteration stops: under RTM the
# per-iteration change ds levels off above the tolerance and the run ends on the no-improvement
# rule; under RDM it reaches the tolerance. Same rung (p=0, T=12), same seed, production settings.
# Reads: data/local/conv_history_sv.jld2  (produced by scripts/analysis/conv_history_sv.jl)
# Run:   julia --project=. scripts/figures/fig_conv.jl
# From a notebook: include this file, then call make_conv() to get the plot object.

ROOT = normpath(joinpath(@__DIR__, "..", ".."))
isdefined(Main, :thesis_plot_theme!) || include(joinpath(ROOT, "src", "thesislib.jl"))
using JLD2, Printf
thesis_plot_theme!()

function make_conv()
    res = load(joinpath(ROOT, "data", "local", "conv_history_sv.jld2"), "res")
    pl = plot(xlabel="iteration", ylabel="\$\\mathrm{d}s\$",
              yscale=:log10, legend=:topright, legendfontsize=8,
              size=thesis_size(0.75; aspect=0.580), margin=2Plots.mm, left_margin=4Plots.mm)
    hline!(pl, [1e-6], color=:black, ls=:dash, lw=0.9, label="tolerance")
    lab(r) = r.reason == "converged" ? "converged at $(r.niters)" : "no improvement, stopped at $(r.niters)"
    for (mode, col, name) in ((:rtm, :crimson, "RTM"), (:rdm, :dodgerblue, "RDM"))
        haskey(res, mode) || continue
        d = res[mode].ds
        plot!(pl, 1:length(d), max.(d, 1e-12), color=col, lw=1.3,
              label="$name, $(lab(res[mode]))")
    end
    for (m, r) in res
        @printf("%s: %s, %d iters, final ds %.2e, min %.2e, lambda0 %.10f%+.10fim\n",
                m, r.reason, r.niters, r.ds[end], minimum(r.ds), real(r.lambda0), imag(r.lambda0))
    end
    return pl
end

if abspath(PROGRAM_FILE) == @__FILE__
    save_thesis_figure(make_conv(), "res_conv")
end

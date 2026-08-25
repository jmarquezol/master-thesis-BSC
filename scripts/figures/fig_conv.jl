# Figure res_conv (thesis fig:conv) — how the block iteration stops: under RTM the change in the
# leading Ritz values levels off above the tolerance and the run ends on the no-improvement rule;
# under RDM it reaches the tolerance.
# Reads: data/local/conv_history.jld2  (produced by scripts/analysis/conv_history.jl)
# Run:   julia --project=. scripts/figures/fig_conv.jl
# From a notebook: include this file, then call make_conv() to get the plot object.

ROOT = normpath(joinpath(@__DIR__, "..", ".."))
isdefined(Main, :thesis_plot_theme!) || include(joinpath(ROOT, "src", "thesislib.jl"))
using JLD2, Printf
thesis_plot_theme!()

function make_conv()
    res = load(joinpath(ROOT, "data", "local", "conv_history.jld2"), "res")
    pl = plot(xlabel="iteration", ylabel="\$\\Delta\\theta\$",
              yscale=:log10, legend=:topright, legendfontsize=8,
              size=thesis_size(0.66; aspect=0.66), margin=2Plots.mm, left_margin=4Plots.mm)
    hline!(pl, [1e-6], color=:black, ls=:dash, lw=0.9, label="tolerance")
    for (mode, col, lab) in ((:rtm, :crimson, "RTM"), (:rdm, :dodgerblue, "RDM"))
        haskey(res, mode) || continue
        d = res[mode].dtheta
        plot!(pl, 1:length(d), max.(d, 1e-12), color=col, lw=1.3,
              label="$lab, $(res[mode].reason) at $(res[mode].niters)")
    end
    for (m, r) in res
        @printf("%s: %s, %d iters, final change %.2e, min %.2e\n",
                m, r.reason, r.niters, r.dtheta[end], minimum(r.dtheta))
    end
    return pl
end

if abspath(PROGRAM_FILE) == @__FILE__
    save_thesis_figure(make_conv(), "res_conv")
end

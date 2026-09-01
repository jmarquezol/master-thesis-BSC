# Figure res_tower (thesis fig:tower) — the boundary operator tower at each coupling: seven
# members read from the phase gaps of the cluster tower arms, against the free-boundary targets.
# Reads: data/cluster/sweep_tower_p{0.0,0.1,0.3,0.5}_bulk.jld2
# Run:   julia --project=. scripts/figures/fig_tower.jl
# From a notebook: include this file, then call make_tower() to get the plot object.

ROOT = normpath(joinpath(@__DIR__, "..", ".."))
isdefined(Main, :thesis_plot_theme!) || include(joinpath(ROOT, "src", "thesislib.jl"))
using JLD2, Printf, Statistics
thesis_plot_theme!()

function make_tower()
    CL = joinpath(ROOT, "data", "cluster")
    V = Dict(0.0 => 2.000, 0.1 => 2.670, 0.3 => 3.967, 0.5 => 5.212)   # v_inf from equilibrium
    arm(file, label) = Dict(k[2] => v for (k, v) in load(joinpath(CL, file), "done")
                            if k[1] == label && !haskey(v, :error))

    # boundary dimensions of one rung: phase gaps to the physical branch, through x = v T dphi/pi
    function dims(e, v)
        i0 = e.i0
        gaps = sort([abs(angle(t / e.theta[i0])) for (i, t) in enumerate(e.theta) if i != i0])
        return v .* e.T .* gaps ./ pi
    end

    sources = [(0.0, "sweep_tower_p0.0_bulk.jld2", "tower_p0.0_bulk"),
               (0.1, "sweep_tower_p0.1_bulk.jld2", "tower_p0.1_bulk"),
               (0.3, "sweep_tower_p0.3_bulk.jld2", "tower_p0.3_bulk"),
               (0.5, "sweep_tower_p0.5_bulk.jld2", "tower_p0.5_bulk")]

    marks = [:circle, :square, :diamond, :utriangle, :dtriangle, :star5, :hexagon]
    cols  = [:dodgerblue, :crimson, :seagreen, :darkorange, :purple, :teal, :brown]
    targets = [0.5, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
    panels = Plots.Plot[]
    for (n, (p, file, label)) in enumerate(sources)
        a = arm(file, label)
        Ts = sort(collect(keys(a)))
        pl = plot(xlabel=(n >= 3 ? "T" : ""), ylabel=(n in (1, 3) ? "x" : ""), legend=false,
                  title="p = $p", ylims=(0.0, 4.8), xlims=(minimum(Ts) - 0.3, maximum(Ts) + 0.3))
        for target in targets
            hline!(pl, [target], color=:black, lw=0.8, ls=:dash, label="")
        end
        for d in 1:7
            ys = [(dd = dims(a[T], V[p]); d <= length(dd) ? dd[d] : NaN) for T in Ts]
            plot!(pl, Ts, ys, color=cols[d], marker=marks[d], ms=3, msw=0, lw=1.2, label="")
        end
        push!(panels, pl)
        @printf("p=%.1f  T=%s\n", p, join([@sprintf("%g", T) for T in Ts], ","))
        @printf("   last rung x1..x7 = %s\n",
                join([@sprintf("%.3f", x) for x in dims(a[Ts[end]], V[p])[1:min(7, end)]], "  "))
    end

    return plot(panels..., layout=(2, 2), size=thesis_size(1.03; aspect=0.551), margin=2Plots.mm)
end

if abspath(PROGRAM_FILE) == @__FILE__
    save_thesis_figure(make_tower(), "res_tower")
end

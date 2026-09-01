# Figure res_bc_pairs (thesis fig:bc_pairs) — the four boundary pairs at p=0: each selects its
# own conformal tower, and the two mixed pairs, which share the same gaps, are told apart by the
# absolute offset x0 read from Eq.(3).
# Reads: data/cluster/sweep_rtm_eigs_p0.0.jld2, data/local/bc_upup_k8.jld2,
#        data/local/bc_updn_k8.jld2, data/local/bc_xup_k8.jld2
# Run:   julia --project=. scripts/figures/fig_bc_pairs.jl
# From a notebook: include this file, then call make_bc_pairs() to get the plot object.

ROOT = normpath(joinpath(@__DIR__, "..", ".."))
isdefined(Main, :thesis_plot_theme!) || include(joinpath(ROOT, "src", "thesislib.jl"))
using JLD2, Printf, Statistics
thesis_plot_theme!()

function make_bc_pairs()
    V0 = 2.0                                   # exact velocity at p=0
    CL = joinpath(ROOT, "data", "cluster")
    LOC = joinpath(ROOT, "data", "local")
    carm(f, l) = Dict(k[2] => v for (k, v) in load(joinpath(CL, f), "done")
                      if k[1] == l && !haskey(v, :error))
    lcache(f) = isfile(f) ? load(f, "done") : Dict{Float64,Any}()
    function byT(d)
        isempty(d) && return Dict{Float64,Any}()
        first(keys(d)) isa Tuple && return Dict{Float64,Any}(k[2] => v for (k, v) in d if k[1] == 0.0)
        return Dict{Float64,Any}(k => v for (k, v) in d)
    end
    # re-anchor the physical branch on modulus continuity in T (never on largest modulus at k=8)
    function reselect(d)
        out = Dict{Float64,Any}(); anchor = nothing
        for T in sort(collect(keys(d)))
            e = d[T]
            i0 = anchor === nothing ? argmax(abs.(e.theta)) : argmin(abs.(abs.(e.theta) .- anchor))
            anchor = abs(e.theta[i0])
            out[T] = (theta=e.theta, i0=i0, theta_phys=e.theta[i0])
        end
        return out
    end

    pairs = [("(a) free--free", reselect(byT(carm("sweep_rtm_eigs_p0.0.jld2", "rtm_eigs_p0.0"))), [0.5, 1.5, 2.0]),
             ("(b) fixed--fixed", reselect(byT(lcache(joinpath(LOC, "bc_upup_k8.jld2")))), [2.0, 3.0, 4.0]),
             ("(c) fixed--antifixed", reselect(byT(lcache(joinpath(LOC, "bc_updn_k8.jld2")))), [1.0, 2.0, 3.0]),
             ("(d) free--fixed", reselect(byT(lcache(joinpath(LOC, "bc_xup_k8.jld2")))), [1.0, 2.0, 3.0])]

    marks = [:circle, :square, :diamond]
    cols  = [:dodgerblue, :crimson, :seagreen]
    panels = Plots.Plot[]
    for (n, (name, d, pred)) in enumerate(pairs)
        Ts = [T for T in sort(collect(keys(d))) if 3.0 <= T <= 9.0]
        pl = plot(xlabel=(n >= 3 ? "T" : ""), ylabel=(n in (1, 3) ? "\$x_i - x_0\$" : ""),
                  legend=false, title=name, titlelocation=:left, titlefontsize=11,
                  ylims=(0.0, maximum(pred) + 0.8), xlims=(2.6, 9.4))
        for target in pred
            hline!(pl, [target], color=:black, lw=0.8, ls=:dash, label="")
        end
        for m in 1:3
            ys = Float64[]
            for T in Ts
                dd = tower_dims(d[T].theta, T, V0; i0=d[T].i0)
                push!(ys, m <= length(dd) ? dd[m] : NaN)
            end
            plot!(pl, Ts, ys, color=cols[m], marker=marks[m], ms=3, msw=0, lw=1.2, label="")
        end
        push!(panels, pl)
        @printf("%-22s rungs T=%g..%g  last = %s\n", name, first(Ts), last(Ts),
                join([@sprintf("%.3f", x) for x in tower_dims(d[last(Ts)].theta, last(Ts), V0; i0=d[last(Ts)].i0)[1:3]], ", "))
    end

    return plot(panels..., layout=(2, 2), size=thesis_size(1.0; aspect=0.570), margin=2Plots.mm)
end

if abspath(PROGRAM_FILE) == @__FILE__
    save_thesis_figure(make_bc_pairs(), "res_bc_pairs")
end

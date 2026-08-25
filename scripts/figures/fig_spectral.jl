# Figures res_circle and res_eq3 (thesis fig:circle, fig:eq3) — the trajectory of mu0 on the
# unit circle, and the Eq.(3) fit per coupling. Arms are merged per coupling as in
# scripts/analysis/cluster_audit.jl.
# Reads: the spectral, tower and entropy arms under data/cluster/ (11 files, see `spec` below)
# Run:   julia --project=. scripts/figures/fig_spectral.jl
# From a notebook: include this file, then call make_spectral(); it returns a named tuple
# (circle=..., eq3=...) with the two plot objects.

ROOT = normpath(joinpath(@__DIR__, "..", ".."))
isdefined(Main, :thesis_plot_theme!) || include(joinpath(ROOT, "src", "thesislib.jl"))
using JLD2, Printf, Statistics, LsqFit
thesis_plot_theme!()

function make_spectral()
    CL = joinpath(ROOT, "data", "cluster")
    V = Dict(0.0 => 2.000, 0.1 => 2.670, 0.3 => 3.967, 0.5 => 5.212)   # v_inf from equilibrium
    arm(file, label) = Dict(k[2] => v for (k, v) in load(joinpath(CL, file), "done")
                            if k[1] == label && !haskey(v, :error))

    spec = Dict{Float64,Dict{Float64,Any}}()
    spec[0.0] = merge(arm("sweep_rtm_eigs_p0.0.jld2", "rtm_eigs_p0.0"),
                      arm("sweep_rtm_eigs_p0.0_fine.jld2", "rtm_eigs_p0.0_fine"),
                      arm("sweep_rtm_eigs_p0.0_fineb.jld2", "rtm_eigs_p0.0_fineb"))
    spec[0.1] = merge(arm("sweep_rtm_eigs_p0.1_bulk.jld2", "rtm_eigs_p0.1_bulk"),
                      arm("sweep_rtm_eigs_p0.1_fine_bulk.jld2", "rtm_eigs_p0.1_fine_bulk"),
                      arm("sweep_tower_p0.1_bulk.jld2", "tower_p0.1_bulk"),
                      arm("sweep_rtm_eigs_p0.1_fineb_bulk.jld2", "rtm_eigs_p0.1_fineb_bulk"))
    spec[0.3] = merge(arm("sweep_ent_p0.3_bulk.jld2", "ent_p0.3_bulk"),
                      arm("sweep_tower_p0.3_bulk.jld2", "tower_p0.3_bulk"),
                      arm("sweep_rtm_eigs_p0.3_fine_bulk.jld2", "rtm_eigs_p0.3_fine_bulk"),
                      arm("sweep_rtm_eigs_p0.3_fineb_bulk.jld2", "rtm_eigs_p0.3_fineb_bulk"))
    # The ladders and windows of Table 1 come from scripts/analysis/extend_window.jl, which selects
    # the rungs beyond the block-method reach. Reading them here keeps the figure and the table
    # showing the same fit.
    SEL = load(joinpath(ROOT, "data", "local", "extended_rungs.jld2"), "selected")
    spec[0.5] = merge(arm("sweep_ent_p0.5_bulk.jld2", "ent_p0.5_bulk"),
                      arm("sweep_tower_p0.5_bulk.jld2", "tower_p0.5_bulk"),
                      arm("sweep_rtm_eigs_p0.5_fine_bulk.jld2", "rtm_eigs_p0.5_fine_bulk"),
                      arm("sweep_rtm_eigs_p0.5_fineb_bulk.jld2", "rtm_eigs_p0.5_fineb_bulk"))
    PS = (0.0, 0.1, 0.3, 0.5)

    # unwrap the phase along the whole ladder before any window is applied
    function unwrapped(mus)
        phi = [angle(-m) for m in mus]
        for i in 2:length(phi)
            d = phi[i] - phi[i-1]
            phi[i] = phi[i-1] + d - 2pi * round(d / (2pi))
        end
        return phi
    end

    # ── figure 1: each coupling on its own circle, drawn at the mean modulus of its ladder, so
    # the different scales stay visible instead of being normalised away
    # the centre of the ring is empty, so the legend sits there rather than stealing width
    pc = plot(xlabel="Re μ₀", ylabel="Im μ₀", legend=:outerright,
              aspect_ratio=:equal, size=thesis_size(0.70; aspect=0.58),
              xlims=(-2.15, 2.15), ylims=(-2.15, 2.15),
              xticks=[-2, -1, 0, 1, 2], yticks=[-2, -1, 0, 1, 2])
    th = range(0, 2pi, length=400)
    for p in PS
        z = SEL[p].mus                      # the same rungs the Eq.(3) fits use
        r = mean(abs.(z))
        plot!(pc, r .* cos.(th), r .* sin.(th), color=P_COLOR[p], ls=:dash, lw=1, alpha=0.6, label="")
        scatter!(pc, real.(z), imag.(z), ms=3, msw=0, color=P_COLOR[p], marker=P_MARKER[p],
                 label=@sprintf("p=%.1f,  |μ₀|=%.3f", p, r))
        @printf("p=%.1f  |mu0| = %.4f  (spread %.2e)\n", p, r, maximum(abs.(z)) - minimum(abs.(z)))
    end

    # ── figure 2: Eq.(3) fit, one panel per coupling, extensive term removed so the curvature shows
    eq3(x, q) = @. q[1] - pi / x + q[2] / x^2
    panels = Plots.Plot[]
    for (n, p) in enumerate(PS)
        Ts, phi = SEL[p].Ts, SEL[p].phis
        y = phi ./ Ts
        f = curve_fit(eq3, Ts, y, [0.1, -0.01])
        c = 24 * V[p] * abs(f.param[2]) / pi
        pl = plot(xlabel=(n >= 3 ? "T" : ""), ylabel=(n in (1, 3) ? "Im λ₀/T − a₀" : ""),
                  legend=:bottomright, title="p = $p")
        scatter!(pl, Ts, y .- f.param[1], ms=3, msw=0, color=P_COLOR[p], label="")
        tt = range(minimum(Ts), maximum(Ts), length=200)
        plot!(pl, tt, eq3(tt, f.param) .- f.param[1], color=:black, lw=1.2,
              label=@sprintf("c = %.2f", c))
        push!(panels, pl)
        @printf("p=%.1f  %2d rungs T=%g..%g  c = %.4f\n", p, length(Ts), first(Ts), last(Ts), c)
    end
    fig = plot(panels..., layout=(2, 2), size=thesis_size(0.95; aspect=0.6), margin=2Plots.mm)

    return (circle=pc, eq3=fig)
end

if abspath(PROGRAM_FILE) == @__FILE__
    sp = make_spectral()
    save_thesis_figure(sp.circle, "res_circle")
    save_thesis_figure(sp.eq3, "res_eq3")
end

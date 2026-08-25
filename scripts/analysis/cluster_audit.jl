ROOT = normpath(joinpath(@__DIR__, "..", ".."))
isdefined(Main, :thesis_plot_theme!) || include(joinpath(ROOT, "src", "thesislib.jl"))
using JLD2, Printf, Statistics, LsqFit

# Audit of the numbers in Table~1 of the thesis.
#
# Sections 1 to 6 recompute the cluster contribution from the cluster caches alone: the block
# iteration's entropy and eigenvalue arms, the dual-unitarity check, and the boundary dimensions.
# p=0 arms are column-independent; every p!=0 arm used here carries the corrected column (_bulk).
#
# The table itself is no longer cluster-only, so Section 7 closes the gap. Its entropy column comes
# from the local single-vector seed ensembles, and its spectral ladders are continued past the
# block iteration's reach by extend_window.jl. Section 7 recomputes both and prints them beside the
# cluster-only values, so the provenance of every quoted number is visible in one place.
const CL = joinpath(ROOT, "data", "cluster")
const V = load(joinpath(ROOT, "data", "local", "alcaraz_velocity.jld2"), "v")

arm(f, l) = Dict(k[2] => v for (k, v) in load(joinpath(CL, f), "done")
                 if k[1] == l && !haskey(v, :error))
ph(z) = angle(-z)
dphw(a, b) = mod(a - b + pi, 2pi) - pi
plateau(s2) = mean(imag.(s2)[max(1, length(s2) ÷ 2 - 1):(length(s2) ÷ 2 + 2)])

ent = Dict(0.0 => merge(arm("sweep_ent_p0.0b_bulk.jld2", "ent_p0.0b_bulk"),
                        arm("sweep_ent_p0.0.jld2", "ent_p0.0")),
           0.1 => merge(arm("sweep_ent_p0.1_bulk.jld2", "ent_p0.1_bulk"),
                        arm("sweep_ent_p0.1b_bulk.jld2", "ent_p0.1b_bulk")),
           0.3 => merge(arm("sweep_ent_p0.3_bulk.jld2", "ent_p0.3_bulk"),
                        arm("sweep_ent_p0.3b_bulk.jld2", "ent_p0.3b_bulk")),
           0.5 => merge(arm("sweep_ent_p0.5_bulk.jld2", "ent_p0.5_bulk"),
                        arm("sweep_ent_p0.5b_bulk.jld2", "ent_p0.5b_bulk")),
           1.0 => arm("sweep_ent_p1.0_bulk_dt0.05.jld2", "ent_p1.0_bulk_dt0.05"))

println("=== 1. entropy plateau, and where each ladder stops ===")
ent_clean = Dict{Float64,Float64}()
for p in (0.0, 0.1, 0.3, 0.5, 1.0)
    Ts = sort(collect(keys(ent[p])))
    pls = [plateau(ent[p][T].s2_base) for T in Ts]
    # the ladder ends where the plateau leaves the family of the ones before it
    lastT = Ts[end]
    for i in 3:length(Ts)
        if abs(pls[i] - median(pls[1:i-1])) > 0.05
            lastT = Ts[i-1]
            break
        end
    end
    ent_clean[p] = lastT
    @printf("p=%.1f  T=%g-%g  clean to %g\n", p, first(Ts), last(Ts), lastT)
    println("   ", join([@sprintf("%g:%.4f", T, x) for (T, x) in zip(Ts, pls)], "  "))
end

println("\n=== 2. the Bou-Comas extrapolation of the plateau, cluster ladders ===")
@. bc(x, q) = q[1] + q[2] / sqrt(x)
@. bc2(x, q) = q[1] + q[2] / x + q[3] / x^2
for p in (0.0, 0.1, 0.3, 0.5)
    Ts = [T for T in sort(collect(keys(ent[p]))) if T <= ent_clean[p]]
    length(Ts) < 4 && (@printf("p=%.1f  only %d rungs, not fitted\n", p, length(Ts)); continue)
    y = [plateau(ent[p][T].s2_base) for T in Ts]
    f1 = curve_fit(bc, Ts, y, [0.1, 0.05])
    @printf("p=%.1f  %2d rungs T=%g-%g   sqrtT: a=%.4f -> c=%.3f",
            p, length(Ts), first(Ts), last(Ts), f1.param[1], 16 * f1.param[1] / pi)
    if length(Ts) >= 5
        f2 = curve_fit(bc2, Ts, y, [0.1, 0.05, 0.0])
        @printf("   two-term: a=%.4f -> c=%.3f", f2.param[1], 16 * f2.param[1] / pi)
    end
    println()
end

println("\n=== 3. Eq.(3) on every corrected-column eigenvalue arm ===")
# spectral arms exist at p=0 and p=0.1; at higher p the tower and entropy arms carry theta_phys
spec = Dict{Float64,Dict{Float64,Any}}()
spec[0.0] = merge(arm("sweep_rtm_eigs_p0.0_fineb.jld2", "rtm_eigs_p0.0_fineb"),
                  arm("sweep_rtm_eigs_p0.0.jld2", "rtm_eigs_p0.0"),
                  arm("sweep_rtm_eigs_p0.0_fine.jld2", "rtm_eigs_p0.0_fine"))
spec[0.1] = merge(arm("sweep_rtm_eigs_p0.1_fineb_bulk.jld2", "rtm_eigs_p0.1_fineb_bulk"),
                  arm("sweep_rtm_eigs_p0.1_bulk.jld2", "rtm_eigs_p0.1_bulk"),
                  arm("sweep_rtm_eigs_p0.1_fine_bulk.jld2", "rtm_eigs_p0.1_fine_bulk"),
                  arm("sweep_tower_p0.1_bulk.jld2", "tower_p0.1_bulk"))
spec[0.3] = merge(ent[0.3], arm("sweep_rtm_eigs_p0.3_fineb_bulk.jld2", "rtm_eigs_p0.3_fineb_bulk"),
                  arm("sweep_tower_p0.3_bulk.jld2", "tower_p0.3_bulk"),
                  arm("sweep_rtm_eigs_p0.3_fine_bulk.jld2", "rtm_eigs_p0.3_fine_bulk"))
spec[0.5] = merge(ent[0.5], arm("sweep_rtm_eigs_p0.5_fineb_bulk.jld2", "rtm_eigs_p0.5_fineb_bulk"),
                  arm("sweep_tower_p0.5_bulk.jld2", "tower_p0.5_bulk"),
                  arm("sweep_rtm_eigs_p0.5_fine_bulk.jld2", "rtm_eigs_p0.5_fine_bulk"))
spec[1.0] = merge(ent[1.0], arm("sweep_tower_p1.0_bulk_dt0.05.jld2", "tower_p1.0_bulk_dt0.05"),
                  arm("sweep_rtm_eigs_p1.0_fine_bulk_dt0.05.jld2", "rtm_eigs_p1.0_fine_bulk_dt0.05"))

function unwrapped(mus)
    p = [ph(m) for m in mus]
    for i in 2:length(p)
        p[i] = p[i-1] + dphw(p[i], p[i-1])
    end
    return p
end
@. eq3(x, q) = q[1] - pi / x + q[2] / x^2

for p in (0.0, 0.1, 0.3, 0.5, 1.0)
    a = spec[p]
    Ts = sort(collect(keys(a)))
    phi = unwrapped([a[T].theta_phys for T in Ts])
    y = phi ./ Ts
    if length(Ts) < 3
        @printf("p=%.1f  %d points only, no fit possible\n", p, length(Ts))
        continue
    end
    f = curve_fit(eq3, Ts, y, [0.1, -0.01])
    c = 24 * V[p] * abs(f.param[2]) / pi
    rms = sqrt(mean((eq3(Ts, f.param) .- y) .^ 2))
    dof = length(Ts) - 2
    @printf("p=%.1f  %2d pts T=%g-%g  dof=%d  |C|=%.5f  c=%.4f  rms=%.1e\n",
            p, length(Ts), first(Ts), last(Ts), dof, abs(f.param[2]), c, rms)
end

println("\n=== 4. modulus spread, the dual-unitarity check ===")
for p in (0.0, 0.1, 0.3, 0.5)
    a = spec[p]
    Ts = sort(collect(keys(a)))
    m = [abs(a[T].theta_phys) for T in Ts]
    @printf("p=%.1f  T=%g-%g  |μ₀| %.4f-%.4f  spread %.2f%%\n",
            p, first(Ts), last(Ts), minimum(m), maximum(m), 100 * (maximum(m) - minimum(m)) / mean(m))
end

println("\n=== 5. x1: free boundary at p=0 from the cluster arm, over its full window ===")
function gap1(e)
    ph0 = ph(e.theta[e.i0])
    best, bg = 0.0, Inf
    for (i, t) in enumerate(e.theta)
        i == e.i0 && continue
        g = dphw(ph(t), ph0)
        abs(g) < bg && (bg = abs(g); best = g)
    end
    return best
end
@. gap_model(T, q) = q[1] / T + q[2] / T^3
free = arm("sweep_rtm_eigs_p0.0.jld2", "rtm_eigs_p0.0")
for hi in (9.0, 14.0, 20.0)
    Ts = [T for T in sort(collect(keys(free))) if T <= hi]
    g = [gap1(free[T]) for T in Ts]
    f = curve_fit(gap_model, Ts, g, [1.0, 0.0])
    xs = V[0.0] .* Ts .* g ./ pi
    @printf("  T=2-%g (%d pts): fit x1=%.4f   per-time %.3f..%.3f  last=%.3f\n",
            hi, length(Ts), V[0.0] * f.param[1] / pi, minimum(xs), maximum(xs), xs[end])
end

println("\n=== 6. x1 at the last time of each corrected-column tower arm ===")
towers = [(0.1, "sweep_tower_p0.1_bulk.jld2", "tower_p0.1_bulk"),
          (0.3, "sweep_tower_p0.3_bulk.jld2", "tower_p0.3_bulk"),
          (0.5, "sweep_tower_p0.5_bulk.jld2", "tower_p0.5_bulk"),
          (1.0, "sweep_tower_p1.0_bulk_dt0.05.jld2", "tower_p1.0_bulk_dt0.05")]
for (p, f, l) in towers
    a = arm(f, l)
    Ts = sort(collect(keys(a)))
    T = last(Ts)
    e = a[T]
    gaps = sort([abs(dphw(ph(t), ph(e.theta[e.i0]))) for (i, t) in enumerate(e.theta) if i != e.i0])
    d = V[p] .* T .* gaps ./ pi
    @printf("p=%.1f  last T=%g  x = %s\n", p, T, join([@sprintf("%.3f", x) for x in d[1:3]], ", "))
end

println("\n=== 7. what Table 1 quotes, and where each column comes from ===")

# The spectral ladders as finally used, written by scripts/analysis/extend_window.jl.
const SELFILE = joinpath(ROOT, "data", "local", "extended_rungs.jld2")
if !isfile(SELFILE)
    println("  extended_rungs.jld2 missing; run scripts/analysis/extend_window.jl first")
else
    sel = load(SELFILE, "selected")
    println("  spectral c: block-iteration rungs continued with single-vector rungs")
    @printf("  %-5s %-12s %-7s %-9s %s\n", "p", "window", "pts", "c", "rms")
    for p in sort(collect(keys(sel)))
        r = sel[p]
        @printf("  %-5.1f T<=%-9g %-7d %-9.4f %.1e\n", p, r.window, length(r.Ts), r.c, r.rms)
    end
end

# The entropy column comes from the local ensembles, not from any cluster arm. The window is the
# range over which the ensembles still return physical plateaus; the surviving fraction is printed
# so the reader can see where each ladder ends and why.
const ENSDIR = joinpath(ROOT, "data", "local", "controls")
const ENT_WINDOW = Dict(0.0 => (8.0, 24.0), 0.1 => (4.0, 13.0),
                        0.3 => (2.0, 9.0), 0.5 => (2.0, 5.0))
physical(x) = 0.05 < x < 0.20
sqrtT(x, q) = @. q[1] + q[2] * x^(-0.5)

println("\n  entropy-route c: local single-vector seed ensembles, median plateau per rung")
@printf("  %-5s %-12s %-7s %-9s %s\n", "p", "window", "rungs", "c", "surviving fraction per rung")
for p in (0.0, 0.1, 0.3, 0.5)
    lo, hi = ENT_WINDOW[p]
    Ts, ys, frac = Float64[], Float64[], String[]
    for T in 2.0:0.5:24.0
        f = joinpath(ENSDIR, "seedens_p$(p)_T$(T).jld2")
        isfile(f) || continue
        ps = [r.plateau for r in values(load(f, "res"))]
        g = filter(physical, ps)
        (lo <= T <= hi) || continue
        isempty(g) && continue
        push!(Ts, T); push!(ys, median(g)); push!(frac, "$(T):$(length(g))/$(length(ps))")
    end
    fit = curve_fit(sqrtT, Ts, ys, [0.098, 0.1])
    @printf("  %-5.1f %-12s %-7d %-9.4f %s\n", p, "$(lo)-$(hi)", length(Ts),
            16 * fit.param[1] / pi, join(frac, " "))
end

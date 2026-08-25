ROOT = normpath(joinpath(@__DIR__, "..", ".."))
isdefined(Main, :thesis_plot_theme!) || include(joinpath(ROOT, "src", "thesislib.jl"))
using JLD2, Printf, Statistics, LsqFit

# Window sensitivity and free-vs-pinned branch constant for the Eq.(3) fits, on the SAME
# round-2 cluster ladders that tab:cp quotes (arm map mirrors cluster_audit.jl).
const CL = joinpath(ROOT, "data", "cluster")
const V = load(joinpath(ROOT, "data", "local", "alcaraz_velocity.jld2"), "v")
arm(f, l) = Dict(k[2] => v for (k, v) in load(joinpath(CL, f), "done")
                 if k[1] == l && !haskey(v, :error))
ph(z) = angle(-z)
dphw(a, b) = mod(a - b + pi, 2pi) - pi

spec = Dict{Float64,Dict{Float64,Any}}()
spec[0.0] = merge(arm("sweep_rtm_eigs_p0.0.jld2", "rtm_eigs_p0.0"),
                  arm("sweep_rtm_eigs_p0.0_fine.jld2", "rtm_eigs_p0.0_fine"))
spec[0.1] = merge(arm("sweep_rtm_eigs_p0.1_bulk.jld2", "rtm_eigs_p0.1_bulk"),
                  arm("sweep_rtm_eigs_p0.1_fine_bulk.jld2", "rtm_eigs_p0.1_fine_bulk"),
                  arm("sweep_tower_p0.1_bulk.jld2", "tower_p0.1_bulk"))
spec[0.3] = merge(arm("sweep_ent_p0.3_bulk.jld2", "ent_p0.3_bulk"),
                  arm("sweep_tower_p0.3_bulk.jld2", "tower_p0.3_bulk"),
                  arm("sweep_rtm_eigs_p0.3_fine_bulk.jld2", "rtm_eigs_p0.3_fine_bulk"))
spec[0.5] = merge(arm("sweep_ent_p0.5_bulk.jld2", "ent_p0.5_bulk"),
                  arm("sweep_tower_p0.5_bulk.jld2", "tower_p0.5_bulk"),
                  arm("sweep_rtm_eigs_p0.5_fine_bulk.jld2", "rtm_eigs_p0.5_fine_bulk"))

function unwrapped(mus)
    p = [ph(m) for m in mus]
    for i in 2:length(p)
        p[i] = p[i-1] + dphw(p[i], p[i-1])
    end
    return p
end
@. eq3_pin(x, q)  = q[1] - pi / x + q[2] / x^2
@. eq3_free(x, q) = q[1] + q[2] / x + q[3] / x^2

for p in (0.0, 0.1, 0.3, 0.5)
    a = spec[p]
    all_Ts = sort(collect(keys(a)))
    # unwrap over the whole ladder first, then window (the branch lives in the full series)
    all_phi = unwrapped([a[T].theta_phys for T in all_Ts])
    println("p=$p  ladder T=$(first(all_Ts))..$(last(all_Ts)) ($(length(all_Ts)) pts)")
    cs = Float64[]
    for Tmin in (2.0, 3.0, 4.0)
        keep = findall(T -> T >= Tmin, all_Ts)
        length(keep) < 5 && continue
        Ts, y = all_Ts[keep], all_phi[keep] ./ all_Ts[keep]
        fp = curve_fit(eq3_pin, Ts, y, [0.1, -0.01])
        ff = curve_fit(eq3_free, Ts, y, [0.1, -pi, -0.01])
        cp = 24V[p] * abs(fp.param[2]) / pi
        cf = 24V[p] * abs(ff.param[3]) / pi
        push!(cs, cp)
        @printf("  T>=%.0f (%2d pts) | pinned c=%.4f rms=%.1e | free B/pi=%+.4f c=%.4f\n",
                Tmin, length(Ts), cp, sqrt(mean(abs2, eq3_pin(Ts, fp.param) .- y)),
                ff.param[2] / pi, cf)
    end
    length(cs) > 1 && @printf("  window sensitivity of pinned c: spread %.3f\n", maximum(cs) - minimum(cs))
end

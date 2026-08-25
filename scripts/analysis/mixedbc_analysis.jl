ROOT = normpath(joinpath(@__DIR__, "..", ".."))
isdefined(Main, :thesis_plot_theme!) || include(joinpath(ROOT, "src", "thesislib.jl"))
using JLD2, Printf, Statistics, LsqFit

# Boundary-pair study at p=0. Four pairs: free-free and fixed-fixed from the cached ladders,
# fixed-antifixed and free-fixed from mixedbc_ladder.jl. Two readings per pair:
#   phase gaps  -> the tower x_i - x_0
#   Eq.(3) fit  -> the absolute offset x_0, which the gaps cannot see
# Branch selection is by continuity in |mu0| along the ladder, never argmax|theta| (at k=8 the
# numerical range of the non-normal operator puts spurious values above the physical one).
const V0 = 2.0
const CL = joinpath(ROOT, "data", "cluster")
carm(f, l) = Dict(k[2] => v for (k, v) in load(joinpath(CL, f), "done")
                  if k[1] == l && !haskey(v, :error))
lcache(f) = isfile(f) ? load(f, "done") : Dict{Float64,Any}()
# bc_upup_k8 is keyed by (p, T); the mixed caches by T
function byT(d)
    isempty(d) && return Dict{Float64,Any}()
    k1 = first(keys(d))
    k1 isa Tuple && return Dict{Float64,Any}(k[2] => v for (k, v) in d if k[1] == 0.0)
    return Dict{Float64,Any}(k => v for (k, v) in d)
end

# reselect i0 by continuity with the previous rung; seed the anchor from the first rung's argmax
function reselect(d)
    Ts = sort(collect(keys(d)))
    out = Dict{Float64,Any}()
    anchor = nothing
    for T in Ts
        e = d[T]
        i0 = anchor === nothing ? argmax(abs.(e.theta)) : argmin(abs.(abs.(e.theta) .- anchor))
        anchor = abs(e.theta[i0])
        out[T] = (theta=e.theta, i0=i0, theta_phys=e.theta[i0])
    end
    return out
end

pairs = [("free-free   (X+,X+)", reselect(byT(carm("sweep_rtm_eigs_p0.0.jld2", "rtm_eigs_p0.0"))), [0.5, 1.5, 2.0], 0.0),
         ("fixed-fixed (Up,Up)", reselect(byT(lcache(joinpath(ROOT, "data", "local", "bc_upup_k8.jld2")))), [2.0, 3.0, 4.0], 0.0),
         ("fixed-anti  (Up,Dn)", reselect(byT(lcache(joinpath(ROOT, "data", "local", "bc_updn_k8.jld2")))), [1.0, 2.0, 3.0], 0.5),
         ("free-fixed  (X+,Up)", reselect(byT(lcache(joinpath(ROOT, "data", "local", "bc_xup_k8.jld2")))), [1.0, 2.0, 3.0], 1/16)]

println("=== 1. the tower of each boundary pair, x_i - x_0 ===")
for (name, d, pred, _) in pairs
    isempty(d) && continue
    Ts = [T for T in sort(collect(keys(d))) if 3.0 <= T <= 9.0]
    isempty(Ts) && continue
    println("\n$name   predicted ", pred)
    for T in Ts
        e = d[T]
        dims = tower_dims(e.theta, T, V0; i0=e.i0)
        @printf("  T=%-4.1f |mu0|=%.6f  x-x0 = %s\n", T, abs(e.theta_phys),
                join([@sprintf("%.3f", x) for x in dims[1:min(3, end)]], ", "))
    end
end

println("\n=== 2. the absolute offset x0 from the leading phase, Eq.(3) ===")
# Im(lam0)/T = a0 + B/T + C/T^2 with C = (pi/v)(x0 - c/24); c = 1/2 exact at p=0
ph(z) = angle(-z)
dphw(a, b) = mod(a - b + pi, 2pi) - pi
function unwrapped(mus)
    p = [ph(m) for m in mus]
    for i in 2:length(p); p[i] = p[i-1] + dphw(p[i], p[i-1]); end
    return p
end
@. eq3_pin(x, q) = q[1] - pi / x + q[2] / x^2
@printf("%-22s %-8s %-10s %-10s %-10s %s\n", "pair", "npts", "C_fit", "x0_fit", "x0_pred", "rms")
for (name, d, _, x0pred) in pairs
    isempty(d) && continue
    Ts = sort(collect(keys(d)))
    length(Ts) < 4 && continue
    phi = unwrapped([d[T].theta_phys for T in Ts])
    y = phi ./ Ts
    f = curve_fit(eq3_pin, Ts, y, [0.1, -0.01])
    C = f.param[2]
    x0 = C * V0 / pi + 0.5 / 24            # C = (pi/v)(x0 - c/24)
    rms = sqrt(mean(abs2, eq3_pin(Ts, f.param) .- y))
    @printf("%-22s %-8d %-10.4f %-10.4f %-10.4f %.1e\n", name, length(Ts), C, x0, x0pred, rms)
end

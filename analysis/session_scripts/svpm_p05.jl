include("../../src/thesislib.jl")
using Printf, JLD2, LsqFit, Statistics

# Single-vector power method on the corrected column at p=0.5, full route to T=10:
# entropy profiles from the dominant pair, lambda0 from the Rayleigh quotient.
p = 0.5
out = "../session_caches/svpm_p05.jld2"
res = isfile(out) ? load(out, "res") : Dict{Float64,Any}()

for T in 2.0:1.0:10.0
    haskey(res, T) && continue
    r = run_pm_diagnosed(T; p=p, maxdim=64, cutoff=1e-12, nbeta=4,
                         itermax=4000, stuck_after=300, column=:bulk5)
    s2 = ITransverse.gen_renyi2(r.L, r.R)
    rig = 1.0 / (norm(r.L) * norm(r.R))
    res[T] = (; lambda0=r.lambda0, s2=collect(s2), rigidity=rig, reason=r.reason, niters=r.niters)
    jldsave(out; res)
    @printf("T=%-4.0f |mu0|=%.6f arg=%.4f  rigidity=%.3g  (%s @%d)\n",
            T, abs(r.lambda0), angle(r.lambda0), rig, r.reason, r.niters)
    flush(stdout)
end

# CFT reads: chord slope and Im plateau per rung (trimmed), |mu0| constancy, Eq.(3) fit
W(t, T) = log((2T / pi) * sin(pi * t / T))
lin(x, q) = q[1] .* x .+ q[2]
trimmed(s2) = collect(s2[3:end-2])
function chord_c(s2, T)
    re = real.(trimmed(s2)); n = length(re)
    ts = range(T/(n+1), T - T/(n+1), length=n); bulk = (n÷4):(3n÷4)
    return 8 * curve_fit(lin, W.(collect(ts[bulk]), T), re[bulk], [0.06, 0.5]).param[1]
end
plateau(s2) = mean(imag.(trimmed(s2))[max(1,end÷2-1):end÷2+2])

println("\nCFT reads per rung (v_inf = 5.212, target c = 1/2, plateau target pi/32 = 0.0982):\n")
@printf("%-4s %-9s %-11s %-9s\n", "T", "c_Re", "Im plateau", "|mu0|")
Ts = sort(collect(keys(res)))
for T in Ts
    @printf("%-4.0f %-9.4f %-11.4f %-9.6f\n", T, chord_c(res[T].s2, T), plateau(res[T].s2), abs(res[T].lambda0))
end

# Eq.(3): Im log(-mu0)/T = a0 + C/T^2, c = 24 v |C| / pi with the free branch constant absorbed
v = 5.21214
lam(T) = log(-res[T].lambda0)
@. eq3(x, q) = q[1] + q[2] / x^2 + q[3] / x
y = [imag(lam(T)) / T for T in Ts]
f = curve_fit(eq3, Ts, y, [0.1, -0.01, -pi])
@printf("\nEq.(3) over T=%g..%g:  C=%.5f -> c = %.3f   (B/pi from 1/T term = %.4f)\n",
        first(Ts), last(Ts), f.param[2], 24 * v * abs(f.param[2]) / pi, f.param[3] / pi)

# Bou-Comas finite-T test on this ladder: does each series decay as T^(-1/2) or 1/T?
@. m_sqrt(x, q) = q[1] + q[2] / sqrt(x)
@. m_lin(x, q)  = q[1] + q[2] / x
println("\nfinite-T decay comparison over the ladder:\n")
for (label, ys) in (("c_Re", [chord_c(res[T].s2, T) for T in Ts]),
                    ("Im plateau", [plateau(res[T].s2) for T in Ts]))
    f1 = curve_fit(m_sqrt, Ts, ys, [0.5, 0.1]); r1 = sqrt(mean(abs2, ys .- m_sqrt(Ts, f1.param)))
    f2 = curve_fit(m_lin,  Ts, ys, [0.5, 0.1]); r2 = sqrt(mean(abs2, ys .- m_lin(Ts, f2.param)))
    conv = label == "Im plateau" ? 16/pi : 1.0
    @printf("%-11s T^-1/2: a=%.4f (c=%.3f) rms=%.2e   1/T: a=%.4f (c=%.3f) rms=%.2e\n",
            label, f1.param[1], conv*f1.param[1], r1, f2.param[1], conv*f2.param[1], r2)
end
println("SVPM-DONE")

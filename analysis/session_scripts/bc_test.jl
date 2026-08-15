using JLD2, LsqFit, Statistics, Printf

# BC finite-T test on the cached corrected-column profiles, with the production convention:
# trim nbeta/2 = 2 cooling bonds per end before any fit.
res = load("../session_caches/entropy_bulkfix.jld2", "results")

W(t, T) = log((2T / pi) * sin(pi * t / T))
lin(x, q) = q[1] .* x .+ q[2]

function trimmed(s2)
    return collect(s2[3:end-2])
end
function chord_c(s2, T)
    re = real.(trimmed(s2)); n = length(re)
    ts = range(T/(n+1), T - T/(n+1), length=n)
    bulk = (n÷4):(3n÷4)
    f = curve_fit(lin, W.(collect(ts[bulk]), T), re[bulk], [0.06, 0.5])
    return 8 * f.param[1]
end
plateau(s2) = mean(imag.(trimmed(s2))[max(1,end÷2-1):end÷2+2])
centre_re(s2) = real(trimmed(s2)[end÷2])

Ts = sort([T for (p, c, T) in keys(res) if p == 0.1 && c == :corrected])
println("corrected column, p=0.1, trimmed convention:\n")
@printf("%-4s %-9s %-11s %-10s\n", "T", "c_Re", "Im plateau", "Re centre")
for T in Ts
    s2 = res[(0.1, :corrected, T)].s2
    @printf("%-4.0f %-9.4f %-11.4f %-10.4f\n", T, chord_c(s2, T), plateau(s2), centre_re(s2))
end

# fit each series with a + b*T^(-1/2) and a + b/T; report extrapolate + rms for both
@. m_sqrt(x, q) = q[1] + q[2] / sqrt(x)
@. m_lin(x, q)  = q[1] + q[2] / x
println("\ndecay-model comparison (4 points, T=2..5):\n")
@printf("%-12s %-22s %-22s\n", "series", "a + b T^-1/2 -> a, rms", "a + b/T -> a, rms")
for (label, f) in (("Im plateau", plateau), ("c_Re", s2 -> 0.0))
    y = label == "Im plateau" ? [plateau(res[(0.1, :corrected, T)].s2) for T in Ts] :
                                [chord_c(res[(0.1, :corrected, T)].s2, T) for T in Ts]
    f1 = curve_fit(m_sqrt, Ts, y, [0.1, 0.1]); r1 = sqrt(mean(abs2, y .- m_sqrt(Ts, f1.param)))
    f2 = curve_fit(m_lin,  Ts, y, [0.1, 0.1]); r2 = sqrt(mean(abs2, y .- m_lin(Ts, f2.param)))
    @printf("%-12s a=%.4f rms=%.2e       a=%.4f rms=%.2e\n", label, f1.param[1], r1, f2.param[1], r2)
end
println("\nIm plateau extrapolations as c = 16a/pi:")
y = [plateau(res[(0.1, :corrected, T)].s2) for T in Ts]
f1 = curve_fit(m_sqrt, Ts, y, [0.1, 0.1]); f2 = curve_fit(m_lin, Ts, y, [0.1, 0.1])
@printf("  T^-1/2 model: c = %.3f    1/T model: c = %.3f\n", 16*f1.param[1]/pi, 16*f2.param[1]/pi)
println("\nproduction column, same fits, for contrast:")
yp = [plateau(res[(0.1, :production, T)].s2) for T in Ts]
g1 = curve_fit(m_sqrt, Ts, yp, [0.1, 0.1]); g2 = curve_fit(m_lin, Ts, yp, [0.1, 0.1])
@printf("  T^-1/2 model: c = %.3f    1/T model: c = %.3f\n", 16*g1.param[1]/pi, 16*g2.param[1]/pi)
println("BC-TEST-DONE")

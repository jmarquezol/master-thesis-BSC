ROOT = normpath(joinpath(@__DIR__, "..", ".."))
isdefined(Main, :thesis_plot_theme!) || include(joinpath(ROOT, "src", "thesislib.jl"))
using JLD2, Printf, Statistics, LsqFit
# Recompute the chi=64 vs chi=128 comparison from the caches, for both estimators.
W(t, T) = log((2T / pi) * sin(pi * t / T))
lin(x, q) = q[1] .* x .+ q[2]
plat(s2) = mean(imag.(s2)[max(1, length(s2) ÷ 2 - 1):(length(s2) ÷ 2 + 2)])
function slope_c(s2, T)
    re = real.(s2); n = length(re)
    ts = range(T / (n + 1), T - T / (n + 1), length=n)
    b = (n ÷ 4):(3n ÷ 4)
    return 8 * curve_fit(lin, W.(collect(ts[b]), T), re[b], [0.06, 0.5]).param[1]
end
trim(s2) = s2[3:end-2]
for (tag, f64, f128) in (("p=0", "svpm_p00.jld2", "svpm_p00_chi128.jld2"),
                         ("p=0.1", "svpm_p01.jld2", "svpm_p01_chi128.jld2"))
    a = load(joinpath(ROOT, "data", "local", f64), "res"); b = load(joinpath(ROOT, "data", "local", f128), "res")
    Ts = sort(intersect(collect(keys(a)), collect(keys(b))))
    isempty(Ts) && continue
    println("\n$tag  shared T: ", Ts)
    @printf("%5s %-10s %-10s %-9s   %-10s %-10s %-9s\n", "T", "plat 64", "plat 128", "d%", "c_Re 64", "c_Re 128", "diff")
    for T in Ts
        p1, p2 = plat(trim(a[T].s2)), plat(trim(b[T].s2))
        c1, c2 = slope_c(trim(a[T].s2), T), slope_c(trim(b[T].s2), T)
        @printf("%5g %-10.4f %-10.4f %-9.2f   %-10.3f %-10.3f %-9.3f\n",
                T, p1, p2, 100abs(p2 - p1) / p1, c1, c2, abs(c2 - c1))
    end
end

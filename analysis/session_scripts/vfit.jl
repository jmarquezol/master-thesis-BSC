using JLD2, Printf, LinearAlgebra
vN = load("../../results/data/nb4_velocity_sizes.jld2", "v")
Ns = [10, 12, 14, 16]
ps = sort(unique(first.(keys(vN))))

println("Finite-size error at p=0, against the exact v=2:\n")
@printf("  %-5s %-10s %-12s %-10s\n", "N", "v(0,N)", "2 - v(0,N)", "N^2 * err")
for N in Ns
    e = 2 - vN[(0.0, N)]
    @printf("  %-5d %-10.5f %-12.5f %-10.4f\n", N, vN[(0.0,N)], e, N^2 * e)
end

println("\nFit v(N) = v_inf + C/N^2 at each p (least squares on the four ring sizes):\n")
@printf("  %-5s %-11s %-9s %-11s %-11s\n", "p", "v_inf", "C", "max resid", "v(N=16)")
for p in ps
    x = [1.0 / N^2 for N in Ns]
    y = [vN[(p, N)] for N in Ns]
    A = hcat(ones(length(x)), x)
    coef = A \ y
    resid = maximum(abs.(y .- A * coef))
    @printf("  %-5.1f %-11.5f %-9.4f %-11.2e %-11.5f\n", p, coef[1], coef[2], resid, vN[(p,16)])
end

println("\nRing size needed for a given accuracy at p=0, from err = C/N^2:")
C0 = 2 - vN[(0.0,16)]
C0 *= 16.0^2
for tol in (0.01, 0.005, 0.001)
    @printf("  |2 - v(0)| < %-6.3f  needs N > %-6.1f  (Hilbert space 2^N = %.3g)\n",
            tol, sqrt(C0/tol), 2.0^ceil(sqrt(C0/tol)))
end

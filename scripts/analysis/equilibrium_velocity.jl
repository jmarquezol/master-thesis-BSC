# Sound velocity from the equilibrium dispersion (Section 4.2, Appendix on the velocity).
#
# For each coupling p we diagonalise the model on periodic rings of N = 10..18 sites, read the
# lowest excitation at zero momentum and at the smallest non-zero momentum, and form
#     v(N) = N [E1(2pi m/N) - E1(0)] / (2 pi |m|),
# the estimator of Alcaraz (2016). The finite-size values are then extrapolated with
# v(N) = v_inf + C/N^2, a form fixed empirically at p = 0 where v = 2 is exact.
#
# Writes: data/local/nb4_velocity_sizes.jld2  (all v(N) per coupling)
#         data/local/alcaraz_velocity.jld2    (the extrapolated v_inf used throughout)
# Run:    julia --project=. scripts/analysis/equilibrium_velocity.jl
#
# This recomputes physics. It is checkpointed per (p, N), so an interrupted run resumes; the
# finished caches ship with the repository and nothing needs to be rerun to reproduce the thesis.

ROOT = normpath(joinpath(@__DIR__, "..", ".."))
isdefined(Main, :thesis_plot_theme!) || include(joinpath(ROOT, "src", "thesislib.jl"))
using JLD2, Printf, LsqFit

const N_LADDER = (10, 12, 14, 16, 18)
const P_VALUES = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0, 1.5]
const SIZES_CACHE = joinpath(ROOT, "data", "local", "nb4_velocity_sizes.jld2")
const VINF_CACHE  = joinpath(ROOT, "data", "local", "alcaraz_velocity.jld2")

# v(N) from one ring, dividing by |m| in case the lowest non-zero momentum is not m = 1
function finite_size_velocity(N, p)
    r = ground_and_gap(N, 1.0, p; howmany=10)
    return N * (r.e1_nonzero - r.e1_zero) / (2 * pi * abs(r.smallest_momentum))
end

vN = isfile(SIZES_CACHE) ? load(SIZES_CACHE, "v") : Dict{Tuple{Float64,Int},Float64}()

for p in P_VALUES, N in N_LADDER
    haskey(vN, (p, N)) && continue
    @printf("ring at p=%.1f, N=%d ...\n", p, N)
    vN[(p, N)] = finite_size_velocity(N, p)
    jldsave(SIZES_CACHE; v=vN)        # save after every ring, not at the end
    GC.gc()
end

# extrapolate: v(N) = v_inf + C/N^2, linear in both parameters
inverse_square(x, q) = @. q[1] + q[2] * x
v_infinity = Dict{Float64,Float64}()

println("\nv(N) = v_inf + C/N^2, fitted over N = 10 to 18\n")
@printf("  %-5s %-9s %-9s %-10s\n", "p", "v_inf", "C", "max resid")
for p in P_VALUES
    x = [1.0 / N^2 for N in N_LADDER]
    y = [vN[(p, N)] for N in N_LADDER]
    fit = curve_fit(inverse_square, x, y, [2.0, -3.0])
    v_infinity[p] = fit.param[1]
    residual = maximum(abs.(y .- inverse_square(x, fit.param)))
    @printf("  %-5.1f %-9.5f %-9.3f %-10.2e\n", p, fit.param[1], fit.param[2], residual)
end

# p = 0 is the control: the exact velocity of the critical Ising chain is 2
@printf("\ncontrol at p=0: v_inf = %.5f against the exact 2, error %.1e\n",
        v_infinity[0.0], abs(2 - v_infinity[0.0]))

jldsave(VINF_CACHE; v=v_infinity)
println("wrote ", relpath(VINF_CACHE, ROOT))

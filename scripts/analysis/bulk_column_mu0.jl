# The test behind the corrected transfer-matrix column of Section 5.
#
# ITransverse builds the column by exponentiating on three sites and repeating the middle tensor.
# With a next-nearest coupling that middle site has no partner two sites away, so the column loses
# a memory channel. We rebuild it on five sites and take the true middle tensor instead
# (column = :bulk5). At p = 0 there is no next-nearest term, so both columns are the same operator.
#
# Showing that the two differ is not enough: it does not say which one is right. So we also compute
# the echo rate of the same model by exact Krylov evolution of an open chain at N = 12, 16, 20,
# extrapolate it linearly in 1/N, and measure both columns against it. The corrected column is the
# one that lands on the exact answer.
#
# Writes: data/local/nb13_mu0.jld2          leading eigenvalue, keyed by (p, column, T)
#         data/local/nb13_krylov_echo.jld2  exact echo rates, keyed by (p, N, T)
# Run:    julia --project=. scripts/analysis/bulk_column_mu0.jl
#
# This recomputes physics; the N = 20 chains are the slow part. Both caches ship with the
# repository and both are checkpointed per point, so an interrupted run resumes.

ROOT = normpath(joinpath(@__DIR__, "..", ".."))
isdefined(Main, :thesis_plot_theme!) || include(joinpath(ROOT, "src", "thesislib.jl"))
using JLD2, Printf, LinearAlgebra, LsqFit
BLAS.set_num_threads(2)

const PS    = (0.0, 0.1, 0.5)
const COLS  = (:legacy3, :bulk5)
const TS    = (1.0, 2.0, 3.0)
const NS    = (12, 16, 20)                    # chain lengths for the exact reference
const MU0_CACHE  = joinpath(ROOT, "data", "local", "nb13_mu0.jld2")
const ECHO_CACHE = joinpath(ROOT, "data", "local", "nb13_krylov_echo.jld2")

# ── 1. the exact reference: echo rate of the open chain, no tensor network involved ──
rates = isfile(ECHO_CACHE) ? load(ECHO_CACHE, "rates") : Dict{Tuple{Float64,Int,Float64},Float64}()
for p in PS, N in NS, T in TS
    haskey(rates, (p, N, T)) && continue
    rates[(p, N, T)] = krylov_echo_rate(N, p, T)
    jldsave(ECHO_CACHE; rates)
    @printf("krylov p=%.1f N=%2d T=%.0f  rate=%.6f\n", p, N, T, rates[(p, N, T)])
end

# the rate still carries a boundary term that falls off as 1/N, so extrapolate it away
linear_in_inverse_N(x, q) = @. q[1] + q[2] * x
exact_rate(p, T) = curve_fit(linear_in_inverse_N, [1 / N for N in NS],
                             [rates[(p, N, T)] for N in NS], [0.4, 0.1]).param[1]

# ── 2. the leading eigenvalue of each column at the same points ──
mu0 = isfile(MU0_CACHE) ? load(MU0_CACHE, "mu0") : Dict{Tuple{Float64,Symbol,Float64},ComplexF64}()
for p in PS, col in COLS, T in TS
    haskey(mu0, (p, col, T)) && continue
    p == 0.0 && col == :bulk5 && continue     # same operator at p=0, the legacy run covers it

    mpo, scaffold = build_alcaraz_tmpo(T; p=p, nbeta=4, column=col)
    theta, _, _, info = block_transfer_eigs(mpo, scaffold; k=2, maxdim=64, cutoff=1e-12,
                                            eigvals_only=true, itermax=300)
    mu0[(p, col, T)] = theta[argmax(abs.(theta))]
    jldsave(MU0_CACHE; mu0)
    @printf("mu0 p=%.1f %-8s T=%.0f  |mu0|=%.6f  (%s)\n", p, col, T, abs(mu0[(p, col, T)]), info[:reason])
end

# ── 3. the comparison: how far each column sits from the exact rate ──
# p=0 is the control. There both columns are the same operator, so their common distance to the
# exact rate is the error floor of the comparison itself; only deviations above it mean anything.
column_rate(p, col, T) = -log(abs(mu0[(p, col, T)]))

println()
@printf("%-5s %-4s %-11s %-13s %-13s\n", "p", "T", "exact", "legacy3 dev", "bulk5 dev")
for p in PS, T in TS
    ell = exact_rate(p, T)
    dev3 = column_rate(p, :legacy3, T) - ell
    dev5 = p == 0.0 ? dev3 : column_rate(p, :bulk5, T) - ell
    @printf("%-5.1f %-4.0f %-11.5f %+-13.5f %+-13.5f\n", p, T, ell, dev3, dev5)
end

println("\nwrote ", relpath(MU0_CACHE, ROOT), " and ", relpath(ECHO_CACHE, ROOT))

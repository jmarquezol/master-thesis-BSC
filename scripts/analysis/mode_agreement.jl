ROOT = normpath(joinpath(@__DIR__, "..", ".."))
isdefined(Main, :thesis_plot_theme!) || include(joinpath(ROOT, "src", "thesislib.jl"))
using JLD2, Printf, Random, LinearAlgebra

# Replaces the unsourced app:blockpm:schur agreement claim: run two INDEPENDENT warm-started
# ladders at p=0.1 (corrected column), one per basis, and compare the Ritz-value ratios.
# :eig  = the de-mixing (eigenvector) mode;  :schur = the eigenvalue-only mode.
BLAS.set_num_threads(2)
const OUT = joinpath(ROOT, "data", "local", "mode_agreement.jld2")
res = isfile(OUT) ? load(OUT, "res") : Dict{Tuple{Symbol,Float64},Any}()

for mode in (:eig, :schur)
    prevL = nothing; prevR = nothing
    for T in 2.0:1.0:8.0
        haskey(res, (mode, T)) && continue
        mpo, scaffold = build_alcaraz_tmpo(T; p=0.1, nbeta=4, column=:bulk5)
        sites = siteinds(scaffold)
        Random.seed!(20260820)
        seedL = prevL === nothing ? nothing : MPS[pad_tmps(v, sites) for v in prevL]
        seedR = prevR === nothing ? nothing : MPS[pad_tmps(v, sites) for v in prevR]
        elapsed = @elapsed theta, L, R, info = block_transfer_eigs(mpo, scaffold;
            k=4, maxdim=64, maxdims=collect(2:2:64),
            cutoff=1e-12, cutoffs=[fill(1e-8, 40); 1e-10],
            itermax=8000, eps_conv=1e-6, itermin=20,
            trunc_mode=:rtm, basis=mode, eigvals_only=(mode === :schur),
            n_track=2, stuck_after=150, seedL=seedL, seedR=seedR)
        i0 = pick_phys_continuity(theta, get(res, (mode, T - 1.0), (theta_phys=nothing,)).theta_phys)
        res[(mode, T)] = (theta=collect(theta), i0=i0, theta_phys=theta[i0],
                          reason=String(info[:reason]), niters=info[:niters], elapsed=elapsed)
        jldsave(OUT; res=res)
        prevL, prevR = L, R
        @printf("%-6s T=%.0f  |mu0|=%.6f  %s %d iters %.0f s\n",
                mode, T, abs(theta[i0]), info[:reason], info[:niters], elapsed)
    end
end

println("\n=== agreement of Ritz-value ratios between the two modes ===")
@printf("%4s  %-12s %-12s %-12s\n", "T", "d|th1/th0|", "d|th2/th0|", "d|th3/th0|")
for T in 2.0:1.0:8.0
    haskey(res, (:eig, T)) && haskey(res, (:schur, T)) || continue
    re_, rs = res[(:eig, T)], res[(:schur, T)]
    se = sort(abs.(re_.theta), rev=true) ./ maximum(abs.(re_.theta))
    ss = sort(abs.(rs.theta), rev=true) ./ maximum(abs.(rs.theta))
    @printf("%4.0f  %-12.2e %-12.2e %-12.2e\n", T, abs(se[2]-ss[2]), abs(se[3]-ss[3]), abs(se[4]-ss[4]))
end

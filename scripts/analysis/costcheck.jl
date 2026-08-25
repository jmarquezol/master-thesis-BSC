ROOT = normpath(joinpath(@__DIR__, "..", ".."))
isdefined(Main, :thesis_plot_theme!) || include(joinpath(ROOT, "src", "thesislib.jl"))
using JLD2, Printf, Random, LinearAlgebra

# The cost subsection claims a speed-up from compressing the partial direct sums, with no
# measurement behind it. `accdim` is the cap on the running sum inside lincomb_mps: at its
# default it equals the production bond dimension, so the sum never grows past chi; setting it
# to k*chi reproduces the plain direct sum, whose peak bond dimension is k times larger.
BLAS.set_num_threads(2)
const OUT = joinpath(ROOT, "data", "local", "costcheck.jld2")
res = isfile(OUT) ? load(OUT, "res") : Dict{Tuple{Symbol,Float64},Any}()

for T in (4.0, 5.0), mode in (:capped, :plain)
    haskey(res, (mode, T)) && continue
    mpo, scaffold = build_alcaraz_tmpo(T; p=0.1, nbeta=4, column=:bulk5)
    Random.seed!(31337 + round(Int, 10T))
    ad = mode === :capped ? 64 : 256          # 64 = chi (production), 256 = k*chi (uncapped)
    elapsed = @elapsed theta, L, R, info = block_transfer_eigs(mpo, scaffold;
        k=4, maxdim=64, maxdims=collect(2:2:64),
        cutoff=1e-12, cutoffs=[fill(1e-8, 40); 1e-10],
        itermax=8000, eps_conv=1e-6, itermin=20,
        trunc_mode=:rtm, basis=:eig, eigvals_only=false,
        n_track=2, stuck_after=150, accdim=ad)
    i0 = argmax(abs.(theta))
    res[(mode, T)] = (theta_phys=theta[i0], reason=String(info[:reason]),
                      niters=info[:niters], elapsed=elapsed, accdim=ad)
    jldsave(OUT; res=res)
    @printf("%-7s T=%.0f accdim=%d  |mu0|=%.8f  %d iters  %.0f s\n",
            mode, T, ad, abs(theta[i0]), info[:niters], elapsed)
end

println("\n=== compressing the partial sums: cost and effect on the answer ===")
for T in (4.0, 5.0)
    haskey(res, (:capped, T)) && haskey(res, (:plain, T)) || continue
    c, p = res[(:capped, T)], res[(:plain, T)]
    @printf("T=%.0f  capped %d its %.0f s | plain %d its %.0f s | ratio %.2f | d|mu0| %.1e relative\n",
            T, c.niters, c.elapsed, p.niters, p.elapsed, p.elapsed / c.elapsed,
            abs(abs(c.theta_phys) - abs(p.theta_phys)) / abs(c.theta_phys))
end

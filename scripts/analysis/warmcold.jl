ROOT = normpath(joinpath(@__DIR__, "..", ".."))
isdefined(Main, :thesis_plot_theme!) || include(joinpath(ROOT, "src", "thesislib.jl"))
using JLD2, Printf, Random, LinearAlgebra

# Warm versus cold start, the one option in the F.4 catalogue with no measurement behind it.
# Same rungs, same settings, same seed; the only difference is whether the block is initialised
# from the converged block of the previous evolution time or from fresh random states.
BLAS.set_num_threads(2)
const OUT = joinpath(ROOT, "data", "local", "warmcold.jld2")
res = isfile(OUT) ? load(OUT, "res") : Dict{Tuple{Symbol,Float64},Any}()
const TS = [3.0, 4.0, 5.0, 6.0]

for mode in (:cold, :warm)
    prevL = nothing; prevR = nothing
    for T in TS
        if haskey(res, (mode, T))
            e = res[(mode, T)]
            continue
        end
        mpo, scaffold = build_alcaraz_tmpo(T; p=0.1, nbeta=4, column=:bulk5)
        sites = siteinds(scaffold)
        Random.seed!(4242 + round(Int, 10T))
        seedL = (mode === :warm && prevL !== nothing) ? MPS[pad_tmps(v, sites) for v in prevL] : nothing
        seedR = (mode === :warm && prevR !== nothing) ? MPS[pad_tmps(v, sites) for v in prevR] : nothing
        elapsed = @elapsed theta, L, R, info = block_transfer_eigs(mpo, scaffold;
            k=4, maxdim=64, maxdims=collect(2:2:64),
            cutoff=1e-12, cutoffs=[fill(1e-8, 40); 1e-10],
            itermax=8000, eps_conv=1e-6, itermin=20,
            trunc_mode=:rtm, basis=:eig, eigvals_only=true,
            n_track=2, stuck_after=150, seedL=seedL, seedR=seedR)
        i0 = argmax(abs.(theta))
        res[(mode, T)] = (theta=collect(theta), theta_phys=theta[i0], reason=String(info[:reason]),
                          niters=info[:niters], elapsed=elapsed)
        jldsave(OUT; res=res)
        if mode === :warm; prevL, prevR = L, R; end
        @printf("%-5s T=%.0f  |mu0|=%.6f  %s  %d iters  %.0f s\n",
                mode, T, abs(theta[i0]), info[:reason], info[:niters], elapsed)
    end
end

println("\n=== warm versus cold at the same rungs ===")
@printf("%5s  %-22s %-22s %s\n", "T", "cold (iters, s)", "warm (iters, s)", "d|mu0| relative")
for T in TS
    haskey(res, (:cold, T)) && haskey(res, (:warm, T)) || continue
    c, w = res[(:cold, T)], res[(:warm, T)]
    d = abs(abs(c.theta_phys) - abs(w.theta_phys)) / abs(c.theta_phys)
    @printf("%5.0f  %-22s %-22s %.2e\n", T,
            @sprintf("%d, %.0f", c.niters, c.elapsed), @sprintf("%d, %.0f", w.niters, w.elapsed), d)
end

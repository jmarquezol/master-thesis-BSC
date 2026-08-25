ROOT = normpath(joinpath(@__DIR__, "..", ".."))
isdefined(Main, :thesis_plot_theme!) || include(joinpath(ROOT, "src", "thesislib.jl"))
using JLD2, Printf, Random, LinearAlgebra

# app:failures states that the failure statistics of the BLOCK iteration (the spectral route) have
# never been measured, only those of the single-vector entropy route. This measures them: the same
# rungs repeated from independent cold seeds, near and past the p=0.3 window edge of T=6.5.
# Lanes write to separate caches so two processes never touch the same file.
BLAS.set_num_threads(2)
const LANE  = parse(Int, ARGS[1])
const SEEDS = LANE == 1 ? (1:3) : (4:6)
const TS    = [6.0, 6.5, 7.0]
const OUT   = joinpath(ROOT, "data", "local", "spectral_seeds_lane$(LANE).jld2")
res = isfile(OUT) ? load(OUT, "res") : Dict{Tuple{Float64,Int},Any}()

for T in TS, s in SEEDS
    haskey(res, (T, s)) && continue
    mpo, scaffold = build_alcaraz_tmpo(T; p=0.3, nbeta=4, column=:bulk5)
    Random.seed!(900_000 + 1000 * s + round(Int, 10T))
    elapsed = @elapsed theta, L, R, info = block_transfer_eigs(mpo, scaffold;
        k=4, maxdim=64, maxdims=collect(2:2:64),
        cutoff=1e-12, cutoffs=[fill(1e-8, 40); 1e-10],
        itermax=8000, eps_conv=1e-6, itermin=20,
        trunc_mode=:rtm, basis=:eig, eigvals_only=true,
        n_track=2, stuck_after=150)
    i0 = argmax(abs.(theta))
    res[(T, s)] = (theta=collect(theta), theta_phys=theta[i0], reason=String(info[:reason]),
                   niters=info[:niters], elapsed=elapsed)
    jldsave(OUT; res=res)
    @printf("T=%.1f seed=%d  |mu0|=%.6f  arg/pi=%+.4f  %s  %d iters  %.0f s\n",
            T, s, abs(theta[i0]), angle(-theta[i0]) / pi, info[:reason], info[:niters], elapsed)
end
println("LANE-$LANE-DONE")

ROOT = normpath(joinpath(@__DIR__, "..", ".."))
isdefined(Main, :thesis_plot_theme!) || include(joinpath(ROOT, "src", "thesislib.jl"))
using JLD2, Printf, Random, LinearAlgebra

# |mu0| agreeing between runs does not mean the GAP agrees: at k=8 two independent fixed-boundary
# ladders once matched on |mu0| to 1e-5 at 13 of 15 rungs while disagreeing on x1 at 7 of them.
# The boundary-pair claim rests on the gaps, so repeat one rung per pair from independent seeds.
BLAS.set_num_threads(2)
const OUT = joinpath(ROOT, "data", "local", "mixedbc_consensus.jld2")
res = isfile(OUT) ? load(OUT, "res") : Dict{Tuple{String,Float64,Int},Any}()

for (pair, bcs) in (("updn", ("Up", "Dn")), ("xup", ("X+", "Up"))), T in (4.0, 5.5), seed in 1:3
    haskey(res, (pair, T, seed)) && continue
    mpo, scaffold = build_tmpo(AlcarazParams(lambda=1.0, p=0.0), AlcarazVD2(), T;
                               dt=0.1, nbeta=4, init_state=bcs[1], init_state_top=bcs[2],
                               column=:legacy3)
    Random.seed!(7000 * seed + round(Int, 10T))
    elapsed = @elapsed theta, L, R, info = block_transfer_eigs(mpo, scaffold;
        k=8, maxdim=64, maxdims=collect(2:2:64),
        cutoff=1e-12, cutoffs=[fill(1e-8, 40); 1e-10],
        itermax=8000, eps_conv=1e-6, itermin=20,
        trunc_mode=:rtm, basis=:eig, eigvals_only=true,
        n_track=2, stuck_after=150)
    res[(pair, T, seed)] = (theta=collect(theta), reason=String(info[:reason]),
                            niters=info[:niters], elapsed=elapsed)
    jldsave(OUT; res=res)
    @printf("%s T=%.1f seed=%d done (%s, %d its, %.0f s)\n", pair, T, seed, info[:reason], info[:niters], elapsed)
end

println("\n=== reproducibility of |mu0| and of the leading gaps ===")
# Anchor the branch on the ladder value at the same rung: at k=8 the numerical range of the
# non-normal operator can place a spurious value above the physical one, so argmax is unsafe.
ladder = Dict("updn" => load(joinpath(ROOT, "data", "local", "bc_updn_k8.jld2"), "done"),
              "xup"  => load(joinpath(ROOT, "data", "local", "bc_xup_k8.jld2"), "done"))
for pair in ("updn", "xup"), T in (4.0, 5.5)
    seeds = [s for s in 1:3 if haskey(res, (pair, T, s))]
    (isempty(seeds) || !haskey(ladder[pair], T)) && continue
    ref = abs(ladder[pair][T].theta_phys)
    @printf("%s T=%.1f  (ladder |mu0| = %.6f)\n", pair, T, ref)
    for s in seeds
        e = res[(pair, T, s)]
        i0 = argmin(abs.(abs.(e.theta) .- ref))
        d = tower_dims(e.theta, T, 2.0; i0=i0)
        @printf("  seed %d  |mu0|=%.6f  x-x0 = %s\n", s, abs(e.theta[i0]),
                join([@sprintf("%.4f", x) for x in d[1:min(3, end)]], ", "))
    end
end

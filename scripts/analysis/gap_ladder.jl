# Qualitative map of the wall: the physical gap |lambda1|/|lambda0| and the peak temporal entropy,
# on a coarse ladder T = 2..10 at p = 0, 0.1, 0.3, 0.5. This is the sweep behind the two panels at
# the start of notebook 03; it shows where the gap closes and the entropy peak grows, not a number
# quoted in the thesis. It uses the three-site column and a small k=4 block, both cheap and both
# adequate for a qualitative picture. Production numbers come from the cluster arms instead.
#
# Writes: data/local/nb55_pgap.jld2
# Run:    julia --project=. scripts/analysis/gap_ladder.jl
#
# This recomputes physics and takes hours. The finished cache ships with the repository. It is
# checkpointed per coupling and warm-started in T, so an interrupted run resumes.

ROOT = normpath(joinpath(@__DIR__, "..", ".."))
isdefined(Main, :thesis_plot_theme!) || include(joinpath(ROOT, "src", "thesislib.jl"))
using JLD2, Printf, LinearAlgebra
BLAS.set_num_threads(2)

const PS    = [0.0, 0.1, 0.3, 0.5]
# p=0 is followed to T=10; the other couplings stop at 7, where the gap has already closed to
# within a couple of per cent and a k=4 block no longer separates the tower reliably.
TS_of(p) = p == 0.0 ? collect(2.0:1.0:10.0) : collect(2.0:1.0:7.0)
const DT    = 0.1
const NBETA = 4
const ALG   = "VD2"
const CACHE = joinpath(ROOT, "data", "local", "nb55_pgap.jld2")

# The spectrum comes in +/- pairs, so lambda0 and its partner are the two largest moduli. We keep
# whichever of the two is closest to the previous step's value, and read the gap off the rest.
function pick_pair(theta, previous)
    previous === nothing && return (1, 2)
    return abs(theta[1] - previous) <= abs(theta[2] - previous) ? (1, 2) : (2, 1)
end

function gap_ratio(theta, i0, ip)
    others = [j for j in eachindex(theta) if j != i0 && j != ip]
    isempty(others) && return 0.0
    return maximum(abs(theta[j]) for j in others) / abs(theta[i0])
end

# One warm-started ladder in T for a single coupling
function gap_ladder(p)
    out = Dict{Float64,Any}()
    prevL = prevR = prev_phys = nothing
    trim = NBETA ÷ 2                       # drop the cooling bonds at each end of the profile

    for T in TS_of(p)
        mpo, scaffold = build_alcaraz_tmpo(T; p=p, lambda=1.0, dt=DT, nbeta=NBETA, MPO_alg=ALG)
        sites = siteinds(scaffold)
        seedL = prevL === nothing ? nothing : MPS[pad_tmps(w, sites) for w in prevL]
        seedR = prevR === nothing ? nothing : MPS[pad_tmps(w, sites) for w in prevR]

        walltime = @elapsed theta, L, R, info = block_transfer_eigs(mpo, scaffold;
            k=4, maxdim=64, maxdims=collect(2:2:64),
            cutoff=1e-12, cutoffs=[fill(1e-8, 40); 1e-10],
            itermax=2000, eps_conv=1e-6, trunc_mode=:rtm,
            n_track=2, stuck_after=150, seedL=seedL, seedR=seedR)

        i0, ip = pick_pair(theta, prev_phys)
        s2 = collect(ITransverse.gen_renyi2(L[i0], R[i0]))[trim+1:end-trim]

        out[T] = (; p, T, dt=DT, MPO_alg=ALG, theta=collect(theta), i0, ip,
                  theta_phys=theta[i0], gap=gap_ratio(theta, i0, ip),
                  s2_re=real.(s2), s2_im=imag.(s2), peak=maximum(real.(s2)),
                  reason=string(info[:reason]), niters=info[:niters],
                  maxchi=max(maxlinkdim(L[i0]), maxlinkdim(R[i0])), walltime)

        @printf("p=%.1f T=%4.1f  gap=%.4f  peak=%.4f  chi=%d  %s in %d iters, %.0f s\n",
                p, T, out[T].gap, out[T].peak, out[T].maxchi, info[:reason], info[:niters], walltime)

        prev_phys = theta[i0]
        prevL, prevR = L, R
        GC.gc()
    end
    return out
end

done = isfile(CACHE) ? load(CACHE, "D") : Dict{Float64,Any}()
for p in PS
    haskey(done, p) && (@printf("p=%.1f cached\n", p); continue)
    done[p] = gap_ladder(p)
    jldsave(CACHE; D=done)                 # checkpoint after every coupling
end
println("wrote ", relpath(CACHE, ROOT))

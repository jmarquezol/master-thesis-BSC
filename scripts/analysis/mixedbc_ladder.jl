ROOT = normpath(joinpath(@__DIR__, "..", ".."))
isdefined(Main, :thesis_plot_theme!) || include(joinpath(ROOT, "src", "thesislib.jl"))
using JLD2, Printf, LinearAlgebra

# Boundary ladders at p=0 for the boundary-pair study of Appendix~app:x1: the fixed-fixed pair
# <Up|..|Up> and the two mixed pairs <Up|..|Dn> and <X+|..|Up>. The fourth pair, free-free, is the
# production cluster arm. All three use the same eigsweep settings so the four are comparable.
#
# Usage:  julia --project=. scripts/analysis/mixedbc_ladder.jl upup|updn|xup
# Writes: data/local/bc_<pair>_k8.jld2
#
# The stored theta_phys is anchored on continuity in T. Do not select the branch by largest
# modulus at k=8: the operator is non-normal, and an unconverged subspace can return a value
# above the physical one. The consumers re-anchor by continuity for the same reason.
# k=8 because the tower gaps are {1,2,3} plus their partners, and a smaller block would cut
# through a degenerate pair. p=0 is the integrable Ising point, so the model has no next-nearest
# term and the three-site column is the correct one here (:bulk5 is needed only for p>0).
BLAS.set_num_threads(2)

const PAIR = ARGS[1]                       # "upup", "updn" or "xup"
const BCS  = Dict("upup" => ("Up", "Up"),
                  "updn" => ("Up", "Dn"),
                  "xup"  => ("X+", "Up"))[PAIR]
# the fixed-fixed tower converges from lower T, so its ladder starts earlier
const TS   = PAIR == "upup" ? collect(2.0:0.5:9.0) : collect(3.0:0.5:8.0)
const CUT  = 1e-8
const KB   = 8
const CACHE = joinpath(ROOT, "data", "local", "bc_$(PAIR)_k$(KB).jld2")
# Keys are the time T. Some shipped caches are keyed by the tuple (p, T); normalise on load so a
# resumed run does not mix the two key types in one file.
function load_done(path)
    isfile(path) || return Dict{Float64,Any}()
    d = load(path, "done")
    isempty(d) && return Dict{Float64,Any}()
    first(keys(d)) isa Tuple && return Dict{Float64,Any}(k[2] => v for (k, v) in d if k[1] == 0.0)
    return Dict{Float64,Any}(k => v for (k, v) in d)
end
done = load_done(CACHE)

previous_L = nothing
previous_R = nothing
for T in TS
    if haskey(done, T)
        @printf("T=%.1f cached\n", T)
        continue
    end
    mpo, scaffold = build_tmpo(AlcarazParams(lambda=1.0, p=0.0), AlcarazVD2(), T;
                               dt=0.1, nbeta=4, init_state=BCS[1], init_state_top=BCS[2],
                               column=:legacy3)
    sites = siteinds(scaffold)
    seedL = previous_L === nothing ? nothing : MPS[pad_tmps(v, sites) for v in previous_L]
    seedR = previous_R === nothing ? nothing : MPS[pad_tmps(v, sites) for v in previous_R]

    elapsed = @elapsed theta, L, R, info = block_transfer_eigs(mpo, scaffold;
        k=KB, maxdim=64, maxdims=collect(2:2:64),
        cutoff=1e-12, cutoffs=[fill(CUT, 40); CUT / 100],
        itermax=8000, eps_conv=1e-6, itermin=20,
        trunc_mode=:rtm, basis=:eig, eigvals_only=true,
        n_track=2, stuck_after=150, seedL=seedL, seedR=seedR)

    i0 = pick_phys_continuity(theta, get(done, T - 0.5, (theta_phys=nothing,)).theta_phys)
    done[T] = (theta=collect(theta), i0=i0, theta_phys=theta[i0],
               reason=String(info[:reason]), niters=info[:niters], elapsed=elapsed)
    jldsave(CACHE; done=done)
    global previous_L, previous_R = L, R
    dims = round.(tower_dims(theta, T, 2.0; i0=i0)[1:min(3, KB - 1)], digits=3)
    @printf("T=%.1f  |mu0|=%.6f  x-x0=%s  %s in %d iters, %.0f s\n",
            T, abs(theta[i0]), string(dims), info[:reason], info[:niters], elapsed)
end

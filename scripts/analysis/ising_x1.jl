# Boundary exponent x1 for the exactly solvable Ising chain, used as the benchmark of notebook 02.
# We run the block power method on the symmetric Ising tMPO for the two boundary conditions and
# read x1 from the first gap of the tower. The exact answers are 1/2 for the free boundary and 2
# for the fixed one, so this checks the whole pipeline against a known result.
#
# A note on the convention: ITransverse builds its Ising model with the field on sigma-z, the
# opposite of the thesis model. In that convention "Up" is the free boundary (x1 = 1/2) and "X+"
# is the fixed one (x1 = 2), which is the reverse of the names used everywhere else here.
#
# Writes: data/local/nb6_ising_x1.jld2
# Run:    julia --project=. scripts/analysis/ising_x1.jl
#
# This recomputes physics and takes about an hour. The finished cache ships with the repository and
# the run is checkpointed per (boundary, T), so an interrupted run resumes.

ROOT = normpath(joinpath(@__DIR__, "..", ".."))
isdefined(Main, :thesis_plot_theme!) || include(joinpath(ROOT, "src", "thesislib.jl"))
using JLD2, Printf, LinearAlgebra
BLAS.set_num_threads(2)

const BCS   = ["Up", "X+"]                 # free, fixed (see the convention note above)
const TS    = [2.0, 3.0, 4.0, 5.0, 6.0, 7.0]
const CACHE = joinpath(ROOT, "data", "local", "nb6_ising_x1.jld2")

# One cold-started block power method run on the Ising tMPO at final time T
function ising_block_eigs(T; init_state, dt=0.1, nbeta=4)
    mpo, scaffold = build_tmpo(IsingParams(1.0, 1.0, 0.0), Murg(), T;
                               dt=dt, nbeta=nbeta, init_state=init_state)
    theta, _, _, info = block_transfer_eigs(mpo, scaffold;
        k=4, maxdim=48, maxdims=collect(2:2:48),
        cutoff=1e-12, cutoffs=[fill(1e-8, 30); 1e-10],
        itermax=2000, eps_conv=1e-6, trunc_mode=:rtm,
        n_track=2, stuck_after=200)
    return collect(theta), info
end

done = isfile(CACHE) ? load(CACHE, "done") : Dict{Tuple{String,Float64},Any}()
for bc in BCS, T in TS
    haskey(done, (bc, T)) && (@printf("%-3s T=%.0f cached\n", bc, T); continue)
    theta, info = ising_block_eigs(T; init_state=bc)
    done[(bc, T)] = (theta=theta, reason=string(info[:reason]), niters=info[:niters])
    jldsave(CACHE; done=done)              # checkpoint after every point
    @printf("%-3s T=%.0f  |theta|=%s  %s in %d iters\n", bc, T,
            join((@sprintf("%.4f", abs(t)) for t in theta), ", "), info[:reason], info[:niters])
    GC.gc()
end
println("wrote ", relpath(CACHE, ROOT))

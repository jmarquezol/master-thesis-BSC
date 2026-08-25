# Cutoff control on the block path: does tightening the singular-value cutoff turn a run that got
# stuck into one that converges, and what does it do to the eigenvalue floor, the gap and the
# plateau? Each point is repeated from three reproducible seeds, at three cutoffs, so the answer
# is an ensemble statement rather than one draw.
#
# The grid covers the two couplings where the block path was hardest (p = 0.3, 0.5) and, as a
# control, the p = 0 fixed boundary, where the answer is known. battery_report.jl reads the caches.
#
# Usage:  julia --project=. scripts/analysis/cutrerun.jl                    # the whole grid
#         julia --project=. scripts/analysis/cutrerun.jl P T CUT SEED [1]   # one point, for lanes
#         (a trailing 1 switches the free boundary |X+> to the fixed |Up>)
# Writes: data/local/controls/cutrerun_p<P>_T<T>_cut<CUT>_s<SEED>[_fixed].jld2
#
# This recomputes physics. The finished caches ship with the repository, and a point that already
# has its cache is skipped, so a rerun costs nothing and an interrupted one resumes.

ROOT = normpath(joinpath(@__DIR__, "..", ".."))
isdefined(Main, :thesis_plot_theme!) || include(joinpath(ROOT, "src", "thesislib.jl"))
using JLD2, Printf, Random, LinearAlgebra
BLAS.set_num_threads(parse(Int, get(ENV, "LANE_BLAS", "2")))

const CACHE_DIR = joinpath(ROOT, "data", "local", "controls")

function cutoff_point(P, T, CUT, SEED; fixed=false)
    bc  = fixed ? "Up" : "X+"
    out = joinpath(CACHE_DIR, "cutrerun_p$(P)_T$(T)_cut$(CUT)_s$(SEED)$(fixed ? "_fixed" : "").jld2")
    isfile(out) && (println("cached: ", basename(out)); return)

    mpo, scaffold = build_tmpo(AlcarazParams(lambda=1.0, p=P), AlcarazVD2(), T;
                               dt=0.1, nbeta=4, init_state=bc, column=:bulk5)
    Random.seed!(round(Int, 1_000_000P) + 1000 * round(Int, 10T) + 100SEED)

    elapsed = @elapsed theta, L, R, info = block_transfer_eigs(mpo, scaffold;
        k=4, maxdim=64, maxdims=collect(2:2:64),
        cutoff=1e-12, cutoffs=[fill(CUT, 40); CUT / 100],
        itermax=8000, eps_conv=1e-6, itermin=20,
        trunc_mode=:rtm, basis=:eig, n_track=2, stuck_after=150)

    i0 = pick_phys_continuity(theta, nothing)
    moduli = sort(abs.(theta), rev=true)
    gap = length(moduli) > 1 ? log(moduli[1] / moduli[2]) : NaN
    s2  = collect(ITransverse.gen_renyi2(L[i0], R[i0]))
    dtheta = info[:dtheta]

    jldsave(out; res=(p=P, T=T, cutoff=CUT, seed=SEED, bc=bc, theta=collect(theta), i0=i0,
                      niters=info[:niters], reason=String(info[:reason]), elapsed=elapsed,
                      g=gap, plateau=plateau_im(s2), s2=s2,
                      bonddim=(maxlinkdim(L[i0]), maxlinkdim(R[i0])),
                      condS=info[:condS], dtheta=collect(dtheta)))

    @printf("p=%.1f T=%.1f cut=%.0e seed=%d %s | %4d its %-9s %6.0fs | g=%.6f plateau=%.4f floor=%.2e\n",
            P, T, CUT, SEED, bc, info[:niters], info[:reason], elapsed, gap,
            plateau_im(s2), dtheta[end])
    flush(stdout)
end

if isempty(ARGS)
    # (p, T, fixed boundary?) -- the two hard couplings, plus the p=0 fixed-boundary control
    points = [(0.3, 3.0, false), (0.3, 4.0, false), (0.5, 3.0, false), (0.0, 3.0, true)]
    for seed in 1:3, cut in (1e-8, 1e-10, 1e-12), (p, T, fixed) in points
        cutoff_point(p, T, cut, seed; fixed=fixed)
    end
else
    cutoff_point(parse(Float64, ARGS[1]), parse(Float64, ARGS[2]),
                 parse(Float64, ARGS[3]), parse(Int, ARGS[4]);
                 fixed=length(ARGS) >= 5 && ARGS[5] == "1")
end
println("cutoff control done")

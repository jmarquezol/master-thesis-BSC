ROOT = normpath(joinpath(@__DIR__, "..", ".."))
isdefined(Main, :thesis_plot_theme!) || include(joinpath(ROOT, "src", "thesislib.jl"))
using LinearAlgebra, Printf, Logging, JLD2, Random
ITensors.disable_warn_order()
Logging.disable_logging(Logging.Warn)

# Exact test of the perturbation argument behind Appendix~app:conditioning, on the corrected
# column. Small enough transfer matrices are built densely and diagonalised outright, so there is
# no power method and no truncation anywhere in the answer. Two questions:
#   1. does the conditioning of the eigenvector basis track the eigenvalue gap, as the appendix says?
#   2. the predicted eigenvalue shift carries a factor 1/rigidity, so at a rigidity of 1e-10 the
#      bound permits an O(1) error -- yet |mu0| is stable to six digits. Is the numerator
#      structurally small, or is the bound simply loose?
# To separate those, each matrix is perturbed twice at equal norm: once by the low-rank error a
# bond-dimension truncation actually makes, and once by an unstructured random matrix as a control.
#
# Usage:  julia --project=. scripts/analysis/dense_bulk5.jl [MAXDIM]
# Writes: data/local/controls/dense_bulk5.jld2
#
# MAXDIM caps the dense matrix size (default 2500); points above it are skipped rather than run.
# The finished cache ships with the repository and each point is saved as it finishes.
BLAS.set_num_threads(parse(Int, get(ENV, "LANE_BLAS", "2")))

const GUARD = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 2500
const OUT   = joinpath(ROOT, "data", "local", "controls", "dense_bulk5.jld2")

function dense_mpo_matrix(mpo::MPO)
    c = mpo[1]
    for i in 2:length(mpo); c = c * mpo[i]; end
    unprimed = [noprime(s) for s in inds(c) if plev(s) == 0]
    col = combiner(unprimed...); row = combiner(prime.(unprimed)...)
    return Matrix(row * c * col, combinedind(row), combinedind(col))
end

# bulk5 column, no cooling sites, so the dense object is the transfer matrix itself
function tmpo_bulk5(p, Nt, dt)
    mpo, scaffold = build_tmpo(AlcarazParams(lambda=1.0, p=p), AlcarazVD2(), Nt * dt;
                               dt=dt, nbeta=0, column=:bulk5)
    s = siteinds(scaffold)
    return mpo, dim(s[1]), length(s), prod(dim.(s))
end

rigid(l, r) = abs(transpose(l) * r) / (norm(l) * norm(r))

# left partner of right eigenvector j, matched by eigenvalue
function pair_at(FR, FL, j)
    iL = argmin(abs.(FL.values .- FR.values[j]))
    return FR.vectors[:, j], FL.vectors[:, iL]
end

res = isfile(OUT) ? load(OUT, "res") : Dict{Any,Any}()

for (p, Nts, dt) in ((0.0, 3:7, 0.5), (0.1, 2:3, 0.5), (0.3, 2:3, 0.5), (0.5, 2:3, 0.5))
    for Nt in Nts
        key = (p, Nt, dt)
        haskey(res, key) && continue
        mpo, d, nsites, dim_dense = tmpo_bulk5(p, Nt, dt)
        dim_dense > GUARD && continue
        @printf("p=%.1f Nt=%d T=%.1f  site dim %d  dense %d ... ", p, Nt, Nt * dt, d, dim_dense)
        flush(stdout)
        t0 = time()
        M  = dense_mpo_matrix(mpo)
        FR = eigen(M); FL = eigen(transpose(M))
        ord = sortperm(abs.(FR.values), rev=true)
        i0, i1 = ord[1], ord[2]
        r0, l0 = pair_at(FR, FL, i0)
        r1, l1 = pair_at(FR, FL, i1)
        gap01  = abs(FR.values[i0] - FR.values[i1])
        dphi   = mod(angle(FR.values[i1]) - angle(FR.values[i0]) + pi, 2pi) - pi

        # perturbation experiment: structured (low-rank truncation of M, i.e. the kind of error a
        # bond-dimension truncation makes) against an unstructured random control at equal norm
        pert = Dict{String,Any}()
        U, S, V = svd(M)
        keepr = max(1, count(>(1e-8 * S[1]), S))
        dE_struct = M - U[:, 1:keepr] * Diagonal(S[1:keepr]) * V[:, 1:keepr]'
        Random.seed!(20260818)
        dE_rand = randn(ComplexF64, size(M)); dE_rand *= norm(dE_struct) / norm(dE_rand)
        for (nm, dE) in (("trunc", dE_struct), ("random", dE_rand))
            pred  = (transpose(l0) * dE * r0) / (transpose(l0) * r0)   # first-order shift
            bound = norm(dE) / rigid(l0, r0)                            # the 1/rigidity bound
            FRp   = eigen(M + dE)
            actual = minimum(abs.(FRp.values .- FR.values[i0]))
            pert[nm] = (; dEnorm=norm(dE), pred=abs(pred), actual=actual, bound=bound,
                          suppression=abs(pred) / bound)
            @printf("\n    %-7s |dE|=%.2e  pred=%.2e  actual=%.2e  bound=%.2e  supp=%.2e",
                    nm, norm(dE), abs(pred), actual, bound, abs(pred) / bound)
        end
        res[key] = (; p, Nt, dt, T=Nt * dt, sitedim=d, nsites, dim=dim_dense,
                      mu0=FR.values[i0], mu1=FR.values[i1], ratio=abs(FR.values[i1] / FR.values[i0]),
                      gap01, dphi_over_pi=dphi / pi,
                      r0=rigid(l0, r0), r1=rigid(l1, r1), condV=cond(FR.vectors),
                      keeprank=keepr, pert, elapsed=time() - t0)
        jldsave(OUT; res)
        @printf("\n    |mu0|=%.5f  |mu1/mu0|=%.5f  dphi/pi=%+.3f  gap=%.3e  r0=%.3e  cond(V)=%.3e  %.0fs\n",
                abs(FR.values[i0]), res[key].ratio, dphi / pi, gap01, res[key].r0,
                res[key].condV, res[key].elapsed)
        flush(stdout)
        M = nothing; FR = nothing; FL = nothing; U = nothing; V = nothing; GC.gc()
    end
end
println("DENSE-DONE")

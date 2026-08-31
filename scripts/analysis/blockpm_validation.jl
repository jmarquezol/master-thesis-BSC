# Validation of the block power method against exact diagonalisation (Appendix app:blockpm:validation).
#
# For each test operator we build the rotated transfer matrix twice: once as an MPO for the block
# iteration, and once as a dense matrix that is diagonalised outright. Comparing the two leading
# moduli isolates the reduced eigensolver from the approximations of the production pipeline.
#
# The three cases answer different questions:
#   1. production column (ANNNI-type, VD2, p = 0.1) — the full pipeline, truncation included;
#   2. a reduced operator small enough that the block never reaches its bond cap, which removes
#      truncation from the comparison and leaves only the reduced eigensolver;
#   3. the critical Ising column, whose transfer matrix is an independently known reference.
#
# Writes data/local/controls/blockpm_validation.jld2 and prints the table of Appendix F.
# Run:  julia --project=. scripts/analysis/blockpm_validation.jl

ROOT = normpath(joinpath(@__DIR__, "..", ".."))
isdefined(Main, :thesis_plot_theme!) || include(joinpath(ROOT, "src", "thesislib.jl"))
using LinearAlgebra, Printf, JLD2, Random, Logging
ITensors.disable_warn_order()
Logging.disable_logging(Logging.Warn)
BLAS.set_num_threads(2)

const OUT = joinpath(ROOT, "data", "local", "controls", "blockpm_validation.jld2")

function dense_mpo_matrix(mpo::MPO)
    c = mpo[1]
    for i in 2:length(mpo); c = c * mpo[i]; end
    unprimed = [noprime(s) for s in inds(c) if plev(s) == 0]
    col = combiner(unprimed...); row = combiner(prime.(unprimed)...)
    return Matrix(row * c * col, combinedind(row), combinedind(col))
end

# one case: build the column, diagonalise it densely, run the block method on the same operator
function validate(label, p, Nt, dt, alg; itermax=300)
    mpo, scaffold = build_tmpo(AlcarazParams(lambda=1.0, p=p), alg, Nt * dt;
                               dt=dt, nbeta=0, column=:bulk5)
    s = siteinds(scaffold)
    sitedim, nsites = dim(s[1]), length(s)
    dense = prod(dim.(s))

    D = dense_mpo_matrix(mpo)
    ev = sort(eigvals(D); by=abs, rev=true)
    exact = abs.(ev[1:2])

    Random.seed!(20260830)
    theta, _, _, info = block_transfer_eigs(mpo, scaffold; k=4, maxdim=64, cutoff=1e-12,
                                            itermax=itermax, eps_conv=1e-10, trunc_mode=:rtm,
                                            basis=:schur, eigvals_only=true)
    got = sort(abs.(theta); rev=true)[1:2]
    err = maximum(abs.(got .- exact))

    @printf("%-34s sites %d x dim %-3d = %-5d   exact |mu0|,|mu1| = %.4f, %.4f   err %.1e  (%s)\n",
            label, nsites, sitedim, dense, exact[1], exact[2], err, info[:reason])
    return (; label, p, Nt, dt, sitedim, nsites, dense, exact, got, err,
            reason=String(info[:reason]))
end

res = isfile(OUT) ? load(OUT, "res") : Dict{String,Any}()

cases = [
    # label                                    p     Nt  dt    construction    itermax
    ("ANNNI-type, p=0.1 (production, VD2)",    0.1,  3,  0.5,  AlcarazVD2(),   300),
    ("reduced operator, p=0 (VD2)",            0.0,  3,  0.5,  AlcarazVD2(),   300),
    ("critical Ising (p=0), 4 sites",          0.0,  4,  0.5,  AlcarazVD2(),  1500),
]

println("=== block method against dense diagonalisation of the same rotated operator ===")
for (label, p, Nt, dt, alg, itmax) in cases
    haskey(res, label) && continue
    res[label] = validate(label, p, Nt, dt, alg; itermax=itmax)
    jldopen(OUT, "w") do f; f["res"] = res; end
end

println("\n--- rows for the appendix table ---")
for (label, _, _, _, _, _) in cases
    haskey(res, label) || continue
    r = res[label]
    @printf("%-38s & \$%d\$ & \$%.4f,\\;%.4f\$ & \$%.1e\$ \\\\\n",
            r.label, r.dense, r.exact[1], r.exact[2], r.err)
end

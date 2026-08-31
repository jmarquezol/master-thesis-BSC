# Producer of the two tables of app:mpo:comparison (tab:mpo_norm, tab:mpo_order).
# A: norm drift of |X+>^N under repeated layers, N=100, p=0.1, dt=0.05, 16 layers.
# B: order of accuracy at fixed physical time T=1, N=12, against exact diagonalisation.
# Run: julia --project=. scripts/analysis/mpo_order_check.jl
ROOT = normpath(joinpath(@__DIR__, "..", ".."))
isdefined(Main, :thesis_plot_theme!) || include(joinpath(ROOT, "src", "thesislib.jl"))
using LinearAlgebra, Printf, Logging
ITensors.disable_warn_order()
Logging.disable_logging(Logging.Warn)

# Two checks on the time-evolution MPO kernels, for app:mpo:comparison.
# A: norm drift of |X+>^N under repeated layers (the existing table, extended to 16 layers).
# B: order of accuracy, from the error at fixed physical time against exact diagonalisation.
BLAS.set_num_threads(4)

function norm_drift(N, p, dt, layers; alg, cutoff=1e-12, maxdim=256)
    sites = siteinds("S=1/2", N)
    U     = expH_alcaraz(sites, 1.0, p; dt=dt, mpo_alg=alg)
    psi   = complex(MPS(sites, "X+"))
    norms = Float64[]
    for _ in 1:layers
        psi = apply(U, psi; cutoff=cutoff, maxdim=maxdim)
        push!(norms, norm(psi))
    end
    return norms, maxlinkdim(U), maxlinkdim(psi)
end

# dense vector of an MPS, in the site ordering of `sites`
function dense_vector(psi::MPS, sites)
    c = psi[1]
    for i in 2:length(psi); c = c * psi[i]; end
    col = combiner(sites...)
    return Vector(Array(c * col, combinedind(col)))
end

function dense_hamiltonian(sites, p)
    H = MPO(alcaraz_opsum(length(sites), 1.0, p), sites)
    c = H[1]
    for i in 2:length(H); c = c * H[i]; end
    col = combiner(sites...); row = combiner(prime.(sites)...)
    return Matrix(row * c * col, combinedind(row), combinedind(col))
end

function order_scan(N, p, T, dts; algs=("WI","WII","VD2"))
    sites = siteinds("S=1/2", N)
    psi0  = complex(MPS(sites, "X+"))
    v0    = dense_vector(psi0, sites)

    Hd = dense_hamiltonian(sites, p)
    F  = eigen(Hermitian(real(Hd)))
    vT = F.vectors * (cis.(-F.values * T) .* (F.vectors' * v0))   # exact e^{-iHT}|psi0>

    errs = Dict(a => Float64[] for a in algs)
    for alg in algs, dt in dts
        nsteps = round(Int, T / dt)
        U   = expH_alcaraz(sites, 1.0, p; dt=dt, mpo_alg=alg)
        psi = deepcopy(psi0)
        for _ in 1:nsteps
            psi = apply(U, psi; cutoff=1e-14, maxdim=1024)
        end
        push!(errs[alg], norm(dense_vector(psi, sites) - vT))
    end
    return errs
end

slope(x, y) = (n = length(x); sx = sum(x); sy = sum(y);
               (n * sum(x .* y) - sx * sy) / (n * sum(abs2, x) - sx^2))

println("=== A. norm drift, N=100, p=0.1, dt=0.05, cutoff 1e-12 ===")
for alg in ("WI", "WII", "VD2")
    norms, dW, dpsi = norm_drift(100, 0.1, 0.05, 16; alg=alg)
    @printf("%-4s  bond(U)=%-3d  maxlinkdim(psi)=%-4d  L1=%.8f  L8=%.8f  L16=%.8f\n",
            alg, dW, dpsi, norms[1], norms[8], norms[16])
end

println("\n=== B. error at fixed T against exact diagonalisation, N=12 ===")
dts = [0.2, 0.1, 0.05, 0.025, 0.0125]
for p in (0.0, 0.1, 0.5)
    errs = order_scan(12, p, 1.0, dts)
    println("p=$p, T=1.0")
    @printf("  %-5s %s\n", "dt", join([@sprintf("%-11s", a) for a in ("WI","WII","VD2")]))
    for (i, dt) in enumerate(dts)
        @printf("  %-5.4f %-11.3e %-11.3e %-11.3e\n", dt, errs["WI"][i], errs["WII"][i], errs["VD2"][i])
    end
    for alg in ("WI", "WII", "VD2")
        sa = slope(log.(dts), log.(errs[alg]))
        st = slope(log.(dts[end-2:end]), log.(errs[alg][end-2:end]))
        @printf("  slope %-4s = %.3f (all)   %.3f (three smallest dt)\n", alg, sa, st)
    end
end

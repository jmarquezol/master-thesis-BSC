# Extension of mpo_order_check.jl part B to smaller time steps: how small must dt be
# before the measured exponent settles on the true order. Source of the p=0.5 numbers and
# of the point where truncation and the reference accuracy limit the measurement.
# Run: julia --project=. scripts/analysis/mpo_dt_tail.jl
ROOT = normpath(joinpath(@__DIR__, "..", ".."))
isdefined(Main, :thesis_plot_theme!) || include(joinpath(ROOT, "src", "thesislib.jl"))
using LinearAlgebra, Printf, Logging
ITensors.disable_warn_order()
Logging.disable_logging(Logging.Warn)

# How small must dt be before the measured error exponent settles on the true order?
# Same set-up as mpo_order_check.jl part B, extended to smaller steps, with the local
# slope between consecutive steps and the wall time of each run.
BLAS.set_num_threads(4)

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

N, T = 12, 1.0
dts  = [0.2, 0.1, 0.05, 0.025, 0.0125, 0.00625, 0.003125]

for p in (0.1, 0.5)
    sites = siteinds("S=1/2", N)
    psi0  = complex(MPS(sites, "X+"))
    v0    = dense_vector(psi0, sites)
    F     = eigen(Hermitian(real(dense_hamiltonian(sites, p))))
    vT    = F.vectors * (cis.(-F.values * T) .* (F.vectors' * v0))

    println("\np=$p, T=$T, N=$N")
    for alg in ("WI", "WII", "VD2")
        errs = Float64[]; secs = Float64[]
        for dt in dts
            t0 = time()
            U   = expH_alcaraz(sites, 1.0, p; dt=dt, mpo_alg=alg)
            psi = deepcopy(psi0)
            for _ in 1:round(Int, T / dt)
                psi = apply(U, psi; cutoff=1e-14, maxdim=1024)
            end
            push!(errs, norm(dense_vector(psi, sites) - vT))
            push!(secs, time() - t0)
        end
        println("  $alg")
        for i in eachindex(dts)
            s = i == 1 ? NaN : log(errs[i-1] / errs[i]) / log(dts[i-1] / dts[i])
            @printf("    dt=%-9.6f err=%-10.3e local slope=%-8s %6.1f s\n",
                    dts[i], errs[i], isnan(s) ? "--" : @sprintf("%.3f", s), secs[i])
        end
    end
end

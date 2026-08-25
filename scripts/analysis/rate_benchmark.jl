# Schrödinger-picture benchmark of the exponential MPO against TDVP: the Loschmidt rate
# l(T) = -log|<psi0|psi(T)>|/N at p=0.1, N=40, from repeated application of U(dt).
# Writes data/local/rate_VD2.jld2 and rate_WII.jld2. The TDVP reference comes from
# tdvp_loschmidt_amplitude, which caches by T in data/local/, so the two are always compared at
# the same times; any time not already cached is evolved on the spot.
# Run:  julia --project=. scripts/analysis/rate_benchmark.jl   (minutes; skips existing caches)

ROOT = normpath(joinpath(@__DIR__, "..", ".."))
isdefined(Main, :thesis_plot_theme!) || include(joinpath(ROOT, "src", "thesislib.jl"))
using JLD2, Printf

const TARGET_TIMES = collect(0.5:1.0:6.5)
const NSITES = 40
const DT = 0.05

function tdvp_rate()
    d = tdvp_loschmidt_amplitude(NSITES, TARGET_TIMES; p=0.1, lambda=1.0, dt=DT)
    return [d[T].rate for T in TARGET_TIMES]
end

function mpo_rate(algorithm)
    cache = joinpath(ROOT, "data", "local", "rate_$(algorithm).jld2")
    isfile(cache) && return load(cache, "rate_$(algorithm)")
    sites = siteinds("S=1/2", NSITES)
    psi0 = complex(MPS(sites, "X+"))
    U = expH_alcaraz(sites, 1.0, 0.1; dt=DT, mpo_alg=algorithm)
    rates = Float64[]
    for T in TARGET_TIMES
        psi = deepcopy(psi0)
        for _ in 1:round(Int, T / DT)
            psi = apply(U, psi; cutoff=1e-14, maxdim=256)
            normalize!(psi)
        end
        push!(rates, -log(abs(inner(psi0, psi))) / NSITES)
        @printf("%s T=%.1f rate=%.5f\n", algorithm, T, rates[end])
    end
    jldopen(cache, "w") do f; f["rate_$(algorithm)"] = rates; end
    return rates
end

if abspath(PROGRAM_FILE) == @__FILE__
    for alg in ("VD2", "WII")
        mpo_rate(alg)
    end
    println("done")
end

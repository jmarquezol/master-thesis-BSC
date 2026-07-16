# cluster/rdm_array/worker.jl
ENV["GKSwstype"] = "100"

include(joinpath(@__DIR__, "..", "..", "src", "thesislib.jl"))

using LinearAlgebra, Printf
using MKL

BLAS.set_num_threads(16)


const P_NNN  = 0.1
const LAMBDA = 1.0
const DT     = 0.1
const NBETA  = 4

function pick_phys(theta_values, previous_physical_value)
    if previous_physical_value === nothing
        return (1, 2)
    end
    distance_to_first  = abs(theta_values[1] - previous_physical_value)
    distance_to_second = abs(theta_values[2] - previous_physical_value)
    if distance_to_first <= distance_to_second
        return (1, 2)
    else
        return (2, 1)
    end
end

function trim_dome(profile, nbeta)
    half = nbeta ÷ 2
    return collect(profile[(half + 1):(end - half)])
end

function phase_rigidity(Lj::MPS, Rj::MPS)
    return 1.0 / (norm(Lj) * norm(Rj))
end

function run_wall_scan_single_T(; T::Float64, chi::Int, label::String,
        cutoff=1e-12, cutoffs=[fill(1e-8, 40); 1e-10],
        trunc_mode=:rtm, basis=:eig,
        itermax=8000, stuck_after=400)
    
    cachefile = joinpath(@__DIR__, "worker_results_T$(T).jld2")
    done = Dict{Tuple{String,Float64},Any}()

    try
        elapsed = @elapsed begin
        mpo, scaffold = build_alcaraz_tmpo(T; p=P_NNN, lambda=LAMBDA, dt=DT, nbeta=NBETA, MPO_alg="VD2")
        
        # No warm starts: seedL and seedR are nothing
        theta, L, R, info = block_transfer_eigs(mpo, scaffold;
            k=4, maxdim=chi, maxdims=collect(2:2:chi),
            cutoff=cutoff, cutoffs=cutoffs,
            itermax=itermax, eps_conv=1e-6, trunc_mode=trunc_mode, basis=basis,
            n_track=2, stuck_after=stuck_after,
            seedL=nothing, seedR=nothing)
        k = length(theta)

        i0, ip = pick_phys(theta, nothing)
        s2_base = trim_dome(ITransverse.gen_renyi2(L[i0], R[i0]), NBETA)

        rigidity = Float64[]
        for j in 1:k
            push!(rigidity, phase_rigidity(L[j], R[j]))
        end
        end # @elapsed

        done[(label, T)] = (label=label, T=T, chi=chi, theta=collect(theta),
            i0=i0, ip=ip, theta_phys=theta[i0], theta_partner=theta[ip],
            s2_base=s2_base, peak=maximum(real.(s2_base)), rigidity=rigidity,
            reason=string(info[:reason]), niters=info[:niters], elapsed=elapsed)

        rigidity_strings = String[]
        for r in rigidity
            push!(rigidity_strings, @sprintf("%.2g", r))
        end
        @info @sprintf("[%s] T=%.1f  %s@%d  |θ0|=%.4f  peak=%.4f  r=[%s]  %.0fs",
            label, T, info[:reason], info[:niters], abs(theta[i0]),
            maximum(real.(s2_base)), join(rigidity_strings, ","), elapsed)
    catch err
        @warn "[$label] T=$T failed: $err"
        done[(label, T)] = (error=string(err),)
    end

    jldsave(cachefile; done=done)
    println("[$label] Saved T=$T to $cachefile")
end

if length(ARGS) < 1
    error("usage: julia worker.jl <task_id>")
end
task_id = parse(Int, ARGS[1])
T_val = 1.0 + task_id

# Hardcoded config for RDM
run_wall_scan_single_T(T=T_val, chi=64, label="rdm64", trunc_mode=:rdm)

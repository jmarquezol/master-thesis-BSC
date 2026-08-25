ROOT = normpath(joinpath(@__DIR__, "..", ".."))
isdefined(Main, :thesis_plot_theme!) || include(joinpath(ROOT, "src", "thesislib.jl"))
using JLD2, Printf, Statistics, LsqFit

# Two parameter verdicts for the F.4 catalogue that need no new simulation, plus the spectral
# seed-ensemble summary. Read-only over existing caches.
const CL = joinpath(ROOT, "data", "cluster")
arm(f, l) = Dict(k[2] => v for (k, v) in load(joinpath(CL, f), "done")
                 if k[1] == l && !haskey(v, :error))
plat(s2) = mean(imag.(s2)[max(1, length(s2) ÷ 2 - 1):(length(s2) ÷ 2 + 2)])

println("=== 1. the cooling regulator beta0: what changes when it moves ===")
# nbeta is encoded in the arm label; beta0 = nbeta*dt/2 at dt=0.1
try
    b = load(joinpath(CL, "sweep_beta_p0.0.jld2"), "done")
    nb_of(k) = parse(Int, match(r"nb(\d+)", k[1]).captures[1])
    Ts = sort(unique([k[2] for k in keys(b)]))
    @printf("%-5s %-6s %-8s %-10s %-10s %-10s\n", "T", "nbeta", "beta0", "|mu0|", "plateau", "corrected")
    for T in Ts[1:min(4, end)]
        ks = sort([k for k in keys(b) if k[2] == T], by=nb_of)
        raw = Float64[]; cor = Float64[]
        for k in ks
            nb = nb_of(k); beta0 = nb * 0.1 / 2
            p0 = plat(b[k].s2_base)
            pc = p0 / (1 - 4beta0 / (pi * T))
            push!(raw, p0); push!(cor, pc)
            @printf("%-5g %-6d %-8.2f %-10.4f %-10.4f %-10.4f\n", T, nb, beta0, abs(b[k].theta_phys), p0, pc)
        end
        sp(x) = 100 * (maximum(x) - minimum(x)) / (sum(x) / length(x))
        @printf("   -> spread across beta0: raw %.1f%%, after dividing by (1-4beta0/piT) %.1f%%\n\n", sp(raw), sp(cor))
    end
catch e
    println("  beta sweep unavailable: ", e)
end

println("\n=== 2. how often the block size had to be escalated ===")
for (p, f, l) in ((0.0, "sweep_tower_p0.0_bulk.jld2", "tower_p0.0_bulk"),
                  (0.1, "sweep_tower_p0.1_bulk.jld2", "tower_p0.1_bulk"),
                  (0.3, "sweep_tower_p0.3_bulk.jld2", "tower_p0.3_bulk"),
                  (0.5, "sweep_tower_p0.5_bulk.jld2", "tower_p0.5_bulk"))
    try
        a = arm(f, l)
        Ts = sort(collect(keys(a)))
        esc = [T for T in Ts if hasproperty(a[T], :escalated) && a[T].escalated]
        ks = unique([a[T].k_used for T in Ts if hasproperty(a[T], :k_used)])
        @printf("p=%.1f  %d rungs, k used %s, escalated at %s\n", p, length(Ts), string(ks),
                isempty(esc) ? "no rung" : join([@sprintf("%g", T) for T in esc], ", "))
    catch e
        @printf("p=%.1f  unavailable\n", p)
    end
end

println("\n=== 3. spectral (block) seed ensemble at p=0.3 ===")
res = Dict{Tuple{Float64,Int},Any}()
for lane in (1, 2)
    f = joinpath(ROOT, "data", "local", "spectral_seeds_lane$(lane).jld2")
    isfile(f) && merge!(res, load(f, "res"))
end
if isempty(res)
    println("  no runs yet")
else
    ref = arm("sweep_rtm_eigs_p0.3_fine_bulk.jld2", "rtm_eigs_p0.3_fine_bulk")
    for T in sort(unique([k[1] for k in keys(res)]))
        seeds = sort([k[2] for k in keys(res) if k[1] == T])
        mus = [abs(res[(T, s)].theta_phys) for s in seeds]
        phs = [angle(-res[(T, s)].theta_phys) / pi for s in seeds]
        refmu = haskey(ref, T) ? abs(ref[T].theta_phys) : NaN
        @printf("T=%-4g %d seeds | |mu0| %.6f..%.6f (spread %.1e) | arg/pi %.4f..%.4f | cluster %.6f\n",
                T, length(seeds), minimum(mus), maximum(mus),
                (maximum(mus) - minimum(mus)) / mean(mus), minimum(phs), maximum(phs), refmu)
    end
end

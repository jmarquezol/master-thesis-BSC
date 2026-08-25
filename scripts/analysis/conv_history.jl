ROOT = normpath(joinpath(@__DIR__, "..", ".."))
isdefined(Main, :thesis_plot_theme!) || include(joinpath(ROOT, "src", "thesislib.jl"))
using JLD2, Printf, Random, LinearAlgebra
# Convergence history for the figure requested in app:failures: the change in the leading Ritz
# values against iteration, under RTM (which levels off above the tolerance) and RDM (which
# reaches it). Same rung, same seed, only the truncation scheme differs.
BLAS.set_num_threads(2)
const OUT = joinpath(ROOT, "data", "local", "conv_history.jld2")
res = isfile(OUT) ? load(OUT, "res") : Dict{Symbol,Any}()
for mode in (:rtm, :rdm)
    haskey(res, mode) && continue
    mpo, scaffold = build_alcaraz_tmpo(4.0; p=0.1, nbeta=4, column=:bulk5)
    Random.seed!(20260822)
    el = @elapsed theta, L, R, info = block_transfer_eigs(mpo, scaffold;
        k=4, maxdim=64, maxdims=collect(2:2:64),
        cutoff=1e-12, cutoffs=[fill(1e-8, 40); 1e-10],
        itermax=8000, eps_conv=1e-6, itermin=20,
        trunc_mode=mode, basis=:eig, eigvals_only=false,
        n_track=2, stuck_after=150)
    i0 = argmax(abs.(theta))
    res[mode] = (dtheta=copy(info[:dtheta]), reason=String(info[:reason]),
                 niters=info[:niters], elapsed=el, theta=collect(theta), theta_phys=theta[i0])
    jldsave(OUT; res=res)
    @printf("%s: %s after %d iters, %.0f s\n", mode, info[:reason], info[:niters], el)
end

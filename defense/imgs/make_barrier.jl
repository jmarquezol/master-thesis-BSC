# Defence-only vector version of the entanglement-barrier plot (slide 4).
#
# Replaces defense/imgs/entanglement_barrier.png, which was a 314x245 raster and
# looked pixelated on a projector. Same data: the maximum bond dimension a TDVP
# run of the p=0.1 chain at N=40 asks for, against evolution time, from
# data/local/tdvp_loschmidt_p0.1_N40.jld2. Nothing in thesis/ is touched.
#
# Run:  julia --project=. defense/imgs/make_barrier.jl     (from the repo root)

ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(ROOT, "src", "thesislib.jl"))
using JLD2, Plots
thesis_plot_theme!()

const CAP = 256

d  = load(joinpath(ROOT, "data", "local", "tdvp_loschmidt_p0.1_N40.jld2"), "done")
Ts = sort(collect(keys(d)))
χs = [d[T].maxchi for T in Ts]

# the run is only meaningful until it reaches the cap
hit = findfirst(==(CAP), χs)

# headroom above the cap so the top-left legend never sits on the dashed line
pl = plot(size = thesis_size(0.58; aspect = 0.78), legend = :topleft,
          xlabel = "evolution time  T", ylabel = "max bond dimension  " * "χ",
          ylims = (0, 400), yticks = 0:100:300,
          xlims = (0.2, 7.3), margin = 2Plots.mm)

hline!(pl, [CAP], color = :gray40, ls = :dash, lw = 1.2, label = "cap χ = $CAP")
plot!(pl, Ts, χs, color = :crimson, lw = 1.6, marker = :square, ms = 3.2, msw = 0,
      label = "TDVP asks for")

out = joinpath(@__DIR__, "entanglement_barrier.pdf")
savefig(pl, out)
println("wrote ", out, "   cap first reached at T = ", Ts[hit])

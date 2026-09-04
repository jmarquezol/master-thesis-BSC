# Defence-only variant of thesis figure res_eq3 (slide 27).
#
# Same data and same fits as scripts/figures/fig_spectral.jl -- this only
# re-sizes the 2x2 panel block so the vertical axis is taller and the curvature
# is legible from the back of a lecture room. The thesis figure is NOT touched;
# \graphicspath in defense/main.tex puts defense/imgs/ first, so this file wins.
#
# Run:  julia --project=. defense/imgs/make_res_eq3_tall.jl     (from the repo root)

ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(ROOT, "scripts", "figures", "fig_spectral.jl"))

const TALL_ASPECT = 0.60          # thesis uses 0.439

sp = make_spectral()
plot!(sp.eq3; size = thesis_size(1.04; aspect = TALL_ASPECT))
out = joinpath(@__DIR__, "res_eq3.pdf")
savefig(sp.eq3, out)
println("wrote ", out, "  size=", thesis_size(1.04; aspect = TALL_ASPECT))

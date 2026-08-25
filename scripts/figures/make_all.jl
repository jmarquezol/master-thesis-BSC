# Regenerate every plot figure of the thesis and sync the PDFs into thesis/imgs/.
# Each figure also lands in figures/ as PNG (quick viewing) and SVG (editable in Inkscape).
# Run:   julia --project=. scripts/figures/make_all.jl

ROOT = normpath(joinpath(@__DIR__, "..", ".."))
isdefined(Main, :thesis_plot_theme!) || include(joinpath(ROOT, "src", "thesislib.jl"))

for f in ("fig_cft_L.jl", "fig_chord_equilibrium.jl", "fig_velocity_vs_p.jl",
          "fig_velocity_extrap.jl", "fig_domes.jl", "fig_domes_hi.jl", "fig_re_c.jl",
          "fig_spectral.jl", "fig_tower.jl", "fig_bc_pairs.jl", "fig_beta0.jl", "fig_conv.jl")
    include(joinpath(@__DIR__, f))
end

# one entry per figure file; fig_spectral.jl produces two figures from one data load
figures = [
    ("cft_L",                 make_cft_L),
    ("cft_chord_equilibrium", make_chord_equilibrium),
    ("velocity_vs_p",         make_velocity_vs_p),
    ("velocity_extrapolation", make_velocity_extrap),
    ("res_domes",             make_domes),
    ("res_domes_hi",          make_domes_hi),
    ("res_re_c",              make_re_c),
    ("res_tower",             make_tower),
    ("res_bc_pairs",          make_bc_pairs),
    ("res_beta0",             make_beta0),
    ("res_conv",              make_conv),
]

names = String[]
for (name, make) in figures
    println("── $name")
    save_thesis_figure(make(), name)
    push!(names, name)
end
println("── res_circle, res_eq3")
sp = make_spectral()
save_thesis_figure(sp.circle, "res_circle")
save_thesis_figure(sp.eq3, "res_eq3")
append!(names, ["res_circle", "res_eq3"])

# the thesis includes the PDF versions
for name in names
    cp(joinpath(ROOT, "figures", "$name.pdf"), joinpath(ROOT, "thesis", "imgs", "$name.pdf");
       force=true)
end
println("synced $(length(names)) PDFs into thesis/imgs/")

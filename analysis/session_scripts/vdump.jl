using JLD2, Printf
vN = load("../../results/data/nb4_velocity_sizes.jld2", "v")
used = load("../../results/data/alcaraz_velocity.jld2", "v")
ps = sort(unique(first.(keys(vN))))
@printf("%-6s %-9s %-9s %-9s %-9s %-11s %-9s\n", "p", "N=10", "N=12", "N=14", "N=16", "in use", "16/12")
for p in ps
    @printf("%-6.1f %-9.4f %-9.4f %-9.4f %-9.4f %-11.4f %-9.5f\n",
            p, vN[(p,10)], vN[(p,12)], vN[(p,14)], vN[(p,16)], used[p], vN[(p,16)]/vN[(p,12)])
end

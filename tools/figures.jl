# Index of the thesis figures: which notebook cell owns each one and which caches it reads.
#
#   julia tools/figures.jl status    what is stale, missing, unreferenced or orphaned
#   julia tools/figures.jl sync      copy the referenced figures into thesis/imgs
#   julia tools/figures.jl regen     re-execute the notebooks owning stale figures, then sync
#
# The plots themselves stay in their owning cell. This file only knows where they come from.

const ROOT     = normpath(joinpath(@__DIR__, ".."))
const NB_DIR   = joinpath(ROOT, "NBs")
const OUT_DIR  = joinpath(ROOT, "results", "imgs")
const TEX_DIR  = joinpath(ROOT, "thesis")
const TEX_IMGS = joinpath(TEX_DIR, "imgs")
const DATA     = joinpath(ROOT, "results", "data")

cache(parts...) = joinpath(DATA, parts...)

const MANIFEST = [
    (img="cft_L.png",                  nb="4_cft_GS.ipynb",           cell="9dd04f64",
     deps=[cache("nb4_fss.jld2")]),
    (img="cft_chord_equilibrium.png",  nb="4_cft_GS.ipynb",           cell="828b1dbc",
     deps=[cache("nb4_chord.jld2")]),
    (img="nb4_velocity_vs_p.png",      nb="4_cft_GS.ipynb",           cell="559f7dd6",
     deps=[cache("alcaraz_velocity.jld2"), cache("nb4_velocity_table.jld2")]),
    (img="gap_closing_wall.png",       nb="5_bPM_gap_closing.ipynb",  cell="c908d781",
     deps=[cache("nb5_rtm_vs_rdm.jld2"), cache("nb3_gap_k6.jld2")]),
    (img="cft_ising_validation.png",   nb="6_loschmidt_ising.ipynb",  cell="a26a4eec",
     deps=[cache("ising_lambda0.jld2"), cache("nb6_ising_x1.jld2")]),
    (img="alcaraz_lambda0_circle.png", nb="7_temp_c.ipynb",           cell="r2_circle",
     deps=[cache("nb8_master.jld2")]),
    (img="alcaraz_tower_spectrum.png", nb="7_temp_c.ipynb",           cell="tower_fig",
     deps=[cache("nb3_p01_stress.jld2")]),
    (img="alcaraz_tower_across_p.png", nb="7_temp_c.ipynb",           cell="tower_acrossp",
     deps=[cache("nb8_master.jld2"), cache("cluster", "sweep_rtm_p0.3.jld2"),
           cache("cluster", "sweep_rtm_p0.5.jld2"),
           cache("cluster", "archive", "sweep_rtm_p0.5_broken.jld2")]),
    (img="eq3_kappa_fit.png",          nb="7_temp_c.ipynb",           cell="eq3_fit",
     deps=[cache("nb8_master.jld2")]),
    (img="val_circle_p0.png",          nb="7_temp_c.ipynb",           cell="p0_spectral_val",
     deps=[cache("nb7_alcaraz_fixedbc.jld2")]),
    (img="val_eq3_p0.png",             nb="7_temp_c.ipynb",           cell="p0_spectral_val",
     deps=[cache("nb7_alcaraz_fixedbc.jld2")]),
    (img="val_x1_p0.png",              nb="7_temp_c.ipynb",           cell="p0_spectral_val",
     deps=[cache("nb7_alcaraz_fixedbc.jld2")]),
    (img="nb13_rigidity.png",          nb="8_eigvec_robustness.ipynb", cell="5d46bafb",
     deps=[cache("nb13_eigvec_ladder.jld2"), cache("nb13_wallscan.jld2")]),
    (img="nb9_psweep_gap.png",         nb="9_cluster_wallscan.ipynb", cell="nb9_s6_fig",
     deps=[cache("cluster", "sweep_rtm_p0.0.jld2"), cache("cluster", "sweep_rtm_p0.1.jld2"),
           cache("cluster", "sweep_rtm_p0.3.jld2")]),
]

# Parsed from the .tex rather than listed, so it cannot drift from what the thesis actually uses.
function referenced_figures()
    names = String[]
    for file in filter(f -> endswith(f, ".tex"), readdir(TEX_DIR))
        for m in eachmatch(r"includegraphics(?:\[[^\]]*\])?\{([^}]*)\}", read(joinpath(TEX_DIR, file), String))
            push!(names, basename(m.captures[1]))
        end
    end
    return sort(unique(names))
end

newest(paths) = maximum([mtime(p) for p in paths if isfile(p)]; init=0.0)

# A cell usually saves its figure just before re-saving its cache, so same-run writes land seconds
# apart in either order. Real drift is hours.
const SAME_RUN = 300.0

function status()
    referenced = referenced_figures()
    println("thesis figures (", length(MANIFEST), " in the manifest, ", length(referenced), " referenced)\n")
    stale = String[]

    for entry in MANIFEST
        built = joinpath(OUT_DIR, entry.img)
        note = if !isfile(built)
            "NOT BUILT"
        elseif mtime(built) < newest(entry.deps) - SAME_RUN
            push!(stale, entry.nb); "STALE (data is newer)"
        elseif !isfile(joinpath(TEX_IMGS, entry.img))
            "not in thesis/imgs"
        elseif mtime(joinpath(TEX_IMGS, entry.img)) < mtime(built) - SAME_RUN
            "thesis/imgs copy is older"
        else
            "ok"
        end
        println("  ", rpad(entry.img, 32), rpad(entry.nb, 26), rpad(entry.cell, 18), note)
    end

    missing_from_manifest = setdiff(referenced, [e.img for e in MANIFEST])
    isempty(missing_from_manifest) || println("\nreferenced but not in the manifest: ", join(missing_from_manifest, ", "))

    unreferenced = setdiff(readdir(TEX_IMGS), referenced)
    isempty(unreferenced) || println("\nunreferenced files in thesis/imgs (", length(unreferenced), "): ",
                                     join(sort(unreferenced), ", "))

    println("\nnotebooks needing a rerun: ", isempty(stale) ? "none" : join(unique(stale), ", "))
    return unique(stale)
end

function sync()
    mkpath(TEX_IMGS)
    for entry in MANIFEST
        built = joinpath(OUT_DIR, entry.img)
        isfile(built) || (println("  skipped (not built): ", entry.img); continue)
        cp(built, joinpath(TEX_IMGS, entry.img); force=true)
        println("  copied ", entry.img)
    end
end

function regen()
    notebooks = status()
    isempty(notebooks) && return println("\nnothing to regenerate")
    for nb in notebooks
        println("\nexecuting ", nb)
        run(`jupyter nbconvert --to notebook --execute --inplace
             --ExecutePreprocessor.timeout=10800 $(joinpath(NB_DIR, nb))`)
    end
    println("\nsyncing")
    sync()
end

mode = isempty(ARGS) ? "status" : ARGS[1]
mode == "status" ? status() :
mode == "sync"   ? sync()   :
mode == "regen"  ? regen()  :
error("usage: julia tools/figures.jl <status|sync|regen>")

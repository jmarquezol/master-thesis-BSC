using JLD2
using Glob

const EXPERIMENT_CACHE = joinpath(@__DIR__, "..", "..", "results", "data", "nb13_wallscan_cutoff.jld2")
const GLOBAL_CACHE = joinpath(@__DIR__, "..", "..", "results", "data", "nb13_wallscan_cluster.jld2")

function collect_results()
    worker_files = glob("worker_results_T*.jld2", @__DIR__)
    if isempty(worker_files)
        println("No worker results found to collect.")
        return
    end

    # Load existing caches if they exist
    exp_done = isfile(EXPERIMENT_CACHE) ? load(EXPERIMENT_CACHE, "done") : Dict{Tuple{String,Float64},Any}()
    global_done = isfile(GLOBAL_CACHE) ? load(GLOBAL_CACHE, "done") : Dict{Tuple{String,Float64},Any}()

    merged_count = 0
    for file in worker_files
        worker_data = load(file, "done")
        for (k, v) in worker_data
            exp_done[k] = v
            global_done[k] = v
            merged_count += 1
        end
    end

    jldsave(EXPERIMENT_CACHE; done=exp_done)
    jldsave(GLOBAL_CACHE; done=global_done)
    println("Successfully collected $merged_count data points.")
    println("Saved to: \n - $EXPERIMENT_CACHE\n - $GLOBAL_CACHE")
end

collect_results()

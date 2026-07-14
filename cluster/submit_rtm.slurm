#!/bin/bash
#SBATCH --chdir .
#SBATCH --job-name=wallscan_rtm-%j
#SBATCH --output=<TODO: absolute path to your cloned repo>/cluster/logs/wallscan_rtm-%j.out
#SBATCH --error=<TODO: absolute path to your cloned repo>/cluster/logs/wallscan_rtm-%j.err
#SBATCH --account=<TODO: your BSC project code, e.g. bsc21>
#SBATCH --qos=<TODO: e.g. gp_bsccase on MN5, ngp_bsccase on GRACE — confirm which system first;
#              use the _debug variant (e.g. gp_debug) for a quick smoke test before the real run>
#SBATCH --time=<TODO: DD-HH:MM:SS, e.g. 2-00:00:00 — T=11..14 is uncharted territory, likely the
#              most expensive points seen yet; request generously, the job is resumable if not>
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=<TODO: e.g. 112 on an MN5 general-purpose node — check what your account
#              can request>

# Re-derived χ=64 baseline (trunc_mode=:rtm), full T=2..14 ladder. Run alongside submit_rdm.slurm
# and submit_cutoff.slurm — all three write into the same results/data/nb13_wallscan_cluster.jld2.

module purge
module load <TODO: BSC's julia module + version, e.g. `module load julia/1.12` — check
#           `module avail julia`>

# This workload is BLAS-bound (ITensors tensor contractions), not Threads.@threads-parallel, so
# point BLAS at the cores SLURM actually gave us rather than relying on Julia's own thread pool.
export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK

cd <TODO: absolute path to your cloned repo>
julia --project=. cluster/wall_scan_cluster.jl rtm

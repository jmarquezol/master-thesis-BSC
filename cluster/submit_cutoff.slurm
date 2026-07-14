#!/bin/bash
#SBATCH --chdir .
#SBATCH --job-name=wallscan_cutoff-%j
#SBATCH --output=<TODO: absolute path to your cloned repo>/cluster/logs/wallscan_cutoff-%j.out
#SBATCH --error=<TODO: absolute path to your cloned repo>/cluster/logs/wallscan_cutoff-%j.err
#SBATCH --account=<TODO: your BSC project code, e.g. bsc21>
#SBATCH --qos=<TODO: e.g. gp_bsccase on MN5, ngp_bsccase on GRACE — confirm which system first;
#              use the _debug variant (e.g. gp_debug) for a quick smoke test before the real run>
#SBATCH --time=<TODO: DD-HH:MM:SS, e.g. 2-00:00:00 — a tighter cutoff keeps more singular values,
#              raising the effective working bond dimension; request generously, resumable>
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=<TODO: e.g. 112 on an MN5 general-purpose node — check what your account
#              can request>

# A tighter cutoff (1e-10/1e-12 instead of the production 1e-8/1e-10) at χ=64, full T=2..14 ladder
# — the other ε knob besides χ itself. Run alongside submit_rtm.slurm and submit_rdm.slurm.

module purge
module load <TODO: BSC's julia module + version, e.g. `module load julia/1.12` — check
#           `module avail julia`>

export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK

cd <TODO: absolute path to your cloned repo>
julia --project=. cluster/wall_scan_cluster.jl cutoff

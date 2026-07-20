#!/bin/bash
#SBATCH --job-name=wallscan_rdm
#SBATCH --output=cluster/logs/wallscan_rdm-%j.out
#SBATCH --error=cluster/logs/wallscan_rdm-%j.err
#SBATCH --account=bsc21
#SBATCH --qos=gp_bsccase
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16

# Warm-started, checkpointed :rdm truncation mode, T=2..12 ladder (the cold array run showed T=9
# alone takes ~20.6h under :rdm; extend Ts in wall_scan_cluster.jl + resubmit if more T is needed
# later — the cache/checkpoint machinery makes that incremental). Run alongside submit_rtm.slurm
# and submit_cutoff.slurm. If this job hits the walltime cap, just `sbatch` it again: it resumes
# from the last checkpointed T, truly warm (not just continuity-anchored).

module load julia/1.12.0

export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK

julia --project=. cluster/wall_scan_cluster.jl rdm

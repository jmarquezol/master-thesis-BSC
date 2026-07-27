#!/bin/bash
#SBATCH --job-name=wallscan_rtm
#SBATCH --output=logs/wallscan_rtm-%j.out
#SBATCH --error=logs/wallscan_rtm-%j.err
#SBATCH --account=bsc21
#SBATCH --qos=gp_bsccase
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16

# Warm-started, checkpointed χ=64 baseline (trunc_mode=:rtm), full T=2..20 ladder (extended from
# T=14 to match the p-sweep arms' reach — this is now the master p=0.1/RTM reference for all of
# them). Run alongside submit_rdm.slurm, submit_cutoff.slurm, and the submit_rtm_p*.slurm p-sweep
# arms — all write into the same results/data/cluster/warm_sweep.jld2, keyed by (label, T), so
# nothing clobbers anything else. If this job hits the walltime cap, just `sbatch` it again: it
# resumes from the last checkpointed T, truly warm (not just continuity-anchored). T=20 is a lot
# further than the original T=14, so expect this job to need at least one resubmission.

# Move to the script's directory so it can be launched from anywhere (matches the *_array/ layout)
cd "${SLURM_SUBMIT_DIR:-.}" || exit 1
mkdir -p logs

module load julia/1.12.0

export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK

julia --project=. cluster/wall_scan_cluster.jl rtm

#!/bin/bash
#SBATCH --job-name=wallscan_cutoff
#SBATCH --output=logs/wallscan_cutoff-%j.out
#SBATCH --error=logs/wallscan_cutoff-%j.err
#SBATCH --account=bsc21
#SBATCH --qos=gp_bsccase
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16

# Warm-started, checkpointed tighter cutoff (1e-10/1e-12 instead of the production 1e-8/1e-10) at
# χ=64, full T=2..14 ladder — the other ε knob besides χ itself. Run alongside submit_rtm.slurm and
# submit_rdm.slurm. If this job hits the walltime cap, just `sbatch` it again: it resumes from the
# last checkpointed T, truly warm (not just continuity-anchored).

cd "${SLURM_SUBMIT_DIR:-.}" || exit 1
mkdir -p logs

module load julia/1.12.0

export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK

julia --project=. cluster/wall_scan_cluster.jl cutoff

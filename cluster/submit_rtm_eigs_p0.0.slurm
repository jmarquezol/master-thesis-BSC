#!/bin/bash
#SBATCH --job-name=eigs_p0.0
#SBATCH --output=logs/eigs_p0.0-%j.out
#SBATCH --error=logs/eigs_p0.0-%j.err
#SBATCH --account=bsc21
#SBATCH --qos=gp_bsccase
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16

# Eigenvalues-only sweep at p=0.0, up to T=20. Only the leading spectrum (no entropy/rigidity),
# so it goes well past the wall. Resubmit this same script if it hits the walltime — it picks up
# from where it left off. Writes to results/data/cluster/sweep_rtm_eigs_p0.0.jld2.

cd "$(dirname "$0")" || exit 1
mkdir -p logs
module load julia/1.12.0
export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK

julia --project=.. wall_scan_cluster.jl eigsweep 0.0 20

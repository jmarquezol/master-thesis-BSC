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

# RTM at p=0.1 with a tighter cutoff, T up to 14. Truncation comparison.

cd "$(dirname "$0")" || exit 1
mkdir -p logs
module load julia/1.12.0
export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
julia --project=.. wall_scan_cluster.jl cutoff

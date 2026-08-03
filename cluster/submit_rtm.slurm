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

# RTM baseline at p=0.1, T up to 20. Part of the original truncation comparison.

cd "$(dirname "$0")" || exit 1
mkdir -p logs
module load julia/1.12.0
export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
julia --project=.. wall_scan_cluster.jl preflight || { echo "PREFLIGHT FAILED - bad environment, not starting the ladder"; exit 1; }

julia --project=.. wall_scan_cluster.jl rtm

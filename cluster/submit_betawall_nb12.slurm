#!/bin/bash
#SBATCH --job-name=betawall_nb12
#SBATCH --output=logs/betawall_nb12-%j.out
#SBATCH --error=logs/betawall_nb12-%j.err
#SBATCH --account=bsc21
#SBATCH --qos=gp_bsccase
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16

# Does the UV regulator move the wall? Test arm: p=0.1 at nbeta=12 (beta0=0.6), T up to 14.
# Compare against submit_betawall_nb4.slurm (beta0=0.2, same code path). Resubmit if it hits
# the walltime -- it resumes warm from the checkpoint.

cd "$(dirname "$0")" || exit 1
mkdir -p logs
module load julia/1.12.0
export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
julia --project=.. wall_scan_cluster.jl preflight || { echo "PREFLIGHT FAILED - bad environment, not starting the ladder"; exit 1; }

julia --project=.. wall_scan_cluster.jl betawall 0.1 12 14

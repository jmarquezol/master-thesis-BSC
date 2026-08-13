#!/bin/bash
#SBATCH --job-name=beta_p0.1
#SBATCH --output=logs/beta_p0.1-%j.out
#SBATCH --error=logs/beta_p0.1-%j.err
#SBATCH --account=bsc21
#SBATCH --qos=gp_bsccase
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16

# β0 regulator scan at p=0.1: the same full run repeated for nbeta = 2, 4, 6, 8
# (β0 = 0.1, 0.2, 0.3, 0.4), up to T=10. Lets us see how much the extracted central charge and
# boundary exponent depend on the amount of imaginary-time cooling. Resubmit if it times out — it
# resumes. Writes to results/data/cluster/sweep_beta_p0.1.jld2.

cd "$(dirname "$0")" || exit 1
mkdir -p logs
module load julia/1.12.0
export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK

julia --project=.. wall_scan_cluster.jl preflight || { echo "PREFLIGHT FAILED - bad environment, not starting the ladder"; exit 1; }

julia --project=.. wall_scan_cluster.jl betascan 0.1 10

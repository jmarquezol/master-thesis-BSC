#!/bin/bash
#SBATCH --job-name=wallscan_rdm
#SBATCH --output=logs/wallscan_rdm-%j.out
#SBATCH --error=logs/wallscan_rdm-%j.err
#SBATCH --account=bsc21
#SBATCH --qos=gp_bsccase
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16

# RDM at p=0.1, T up to 12 (about where the RTM p=0.1 run reached). One job: does switching the
# truncation from RTM to RDM move the wall at a frustrated point? Compare vs sweep_rtm_p0.1.jld2.
# RDM is a few times slower than RTM. Writes to results/data/cluster/sweep_rdm_p0.1.jld2.

cd "$(dirname "$0")" || exit 1
mkdir -p logs
module load julia/1.12.0
export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
julia --project=.. wall_scan_cluster.jl rdm

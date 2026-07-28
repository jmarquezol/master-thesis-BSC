#!/bin/bash
#SBATCH --job-name=wallscan_rtm_p0.5
#SBATCH --output=logs/wallscan_rtm_p0.5-%j.out
#SBATCH --error=logs/wallscan_rtm_p0.5-%j.err
#SBATCH --account=bsc21
#SBATCH --qos=gp_bsccase
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16

# Full run (with entropy) at p=0.5, RTM, T up to 20. Resubmit if it times out — it resumes.

cd "${SLURM_SUBMIT_DIR:-.}" || exit 1
mkdir -p logs
module load julia/1.12.0
export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
julia --project=.. wall_scan_cluster.jl psweep 0.5 20

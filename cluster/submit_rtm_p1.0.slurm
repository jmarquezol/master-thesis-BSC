#!/bin/bash
#SBATCH --job-name=wallscan_rtm_p1.0
#SBATCH --output=logs/wallscan_rtm_p1.0-%j.out
#SBATCH --error=logs/wallscan_rtm_p1.0-%j.err
#SBATCH --account=bsc21
#SBATCH --qos=gp_bsccase
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16

# p-sweep arm: p=1.0, RTM route, warm-started and checkpointed, all the way to T=20. See
# submit_rtm_p0.0.slurm for the full explanation — same design, just a different NNN coupling.

cd "$(dirname "$0")" || exit 1
mkdir -p logs

module load julia/1.12.0

export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK

julia --project=.. wall_scan_cluster.jl psweep 1.0 20

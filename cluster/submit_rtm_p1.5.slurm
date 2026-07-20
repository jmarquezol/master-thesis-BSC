#!/bin/bash
#SBATCH --job-name=wallscan_rtm_p1.5
#SBATCH --output=cluster/logs/wallscan_rtm_p1.5-%j.out
#SBATCH --error=cluster/logs/wallscan_rtm_p1.5-%j.err
#SBATCH --account=bsc21
#SBATCH --qos=gp_bsccase
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16

# p-sweep arm: p=1.5, RTM route, warm-started and checkpointed, all the way to T=20. See
# submit_rtm_p0.0.slurm for the full explanation — same design, just a different NNN coupling.
# p=1.5 is the largest NNN coupling in this campaign (CLAUDE.md's DMRG sweep covers p<=~2 as the
# confirmed Ising-class window) — the one most likely to need the k=4->6 escalation.


module load julia/1.12.0

export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK

julia --project=. cluster/wall_scan_cluster.jl psweep 1.5 20

#!/bin/bash
#SBATCH --job-name=wallscan_rtm_p0.0
#SBATCH --output=cluster/logs/wallscan_rtm_p0.0-%j.out
#SBATCH --error=cluster/logs/wallscan_rtm_p0.0-%j.err
#SBATCH --account=bsc21
#SBATCH --qos=gp_bsccase
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40

# v3 p-sweep arm: p=0.0, RTM route, TIGHT cutoff (1e-12/1e-14), k=4 (no k=6 escalation),
# warm-started + checkpointed, capped at T=14. Writes its OWN cache
# results/data/cluster/sweep_rtm_p0.0.jld2 (per-label -> no more shared-file clobbering, the bug
# that lost most of the v2 compute). Walltime cap loses no progress: just resubmit, it resumes from
# the last checkpointed T. --mem=120G fixes the v2 OOM kills (raise toward --mem=0 = whole node if
# OOM recurs). When all arms are done, run: julia --project=. cluster/merge_sweeps.jl

module load julia/1.12.0

export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK

julia --project=. cluster/wall_scan_cluster.jl psweep 0.0 14

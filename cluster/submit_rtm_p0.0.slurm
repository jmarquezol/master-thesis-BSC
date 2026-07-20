#!/bin/bash
#SBATCH --job-name=wallscan_rtm_p0.0
#SBATCH --output=cluster/logs/wallscan_rtm_p0.0-%j.out
#SBATCH --error=cluster/logs/wallscan_rtm_p0.0-%j.err
#SBATCH --account=bsc21
#SBATCH --qos=gp_bsccase
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16

# p-sweep arm: p=0.0 (pure TFIM, the Ising calibration point), RTM route, warm-started and
# checkpointed, all the way to T=20 — much further than the original p=0.1 triple (T=14), so this
# job is very likely to need several resubmissions before it finishes. That is expected and safe:
# hitting the walltime cap loses no progress, just `sbatch` this same script again and it resumes
# from the last checkpointed T, truly warm (not just continuity-anchored). Runs alongside the other
# psweep arms (p=0.3/0.5/1.0/1.5) and the original p=0.1 triple — all write into the same
# results/data/cluster/warm_sweep.jld2, keyed by (label, T), so nothing clobbers anything else.


module load julia/1.12.0

export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK

julia --project=. cluster/wall_scan_cluster.jl psweep 0.0 20

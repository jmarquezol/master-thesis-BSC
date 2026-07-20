#!/bin/bash
#SBATCH --job-name=wallscan_rdm
#SBATCH --output=cluster/logs/wallscan_rdm-%j.out
#SBATCH --error=cluster/logs/wallscan_rdm-%j.err
#SBATCH --account=bsc21
#SBATCH --qos=gp_bsccase
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G

# v3 RDM diagnostic — RUN SEPARATELY, and only AFTER the 5 RTM arms have landed clean. It is a
# cross-check, not a result: NB9 already found RDM buys no physics over RTM for 4-11x the cost.
# Own cache results/data/cluster/sweep_rdm64.jld2 (per-label), --mem=120G, T=12.
#
# IMPORTANT: v2's rdm job produced ZERO output in 2 full days (only the shell module-load line) and
# was walltime-killed. Likely causes now addressed: the shared-file stall (fixed by per-label
# caches) or the RDM path hanging. WATCH the first ~10 minutes of cluster/logs/wallscan_rdm-*.err:
# if it prints T=2 progress, let it run; if it is still silent after ~10 min, KILL it and report
# rather than burning another 2 days.

module load julia/1.12.0

export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK

julia --project=. cluster/wall_scan_cluster.jl rdm

#!/bin/bash
#SBATCH --job-name=wallscan_rtm_p1.0
#SBATCH --output=cluster/logs/wallscan_rtm_p1.0-%j.out
#SBATCH --error=cluster/logs/wallscan_rtm_p1.0-%j.err
#SBATCH --account=bsc21
#SBATCH --qos=gp_bsccase
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40

# v3 REACH arm: p=1.0 (strong frustration). The tight cutoff (1e-12/1e-14) prevents the v2
# null-collapse (v2 gave |θ0|=0 here; v3 gives finite, growing |θ0|). CAVEAT: |θ0| is large and
# grows with p (a non-universal normalisation effect) and the spectrum's top is near-degenerate
# (tiny gaps) — so λ0-selection is fragile and these arms are EXPLORATORY. The phase-based CFT
# extractions (Eq.3 velocity, Eq.4 x1) are modulus-independent and may still be valid; assess in
# post-processing. RTM route, k=4, warm-started + checkpointed, T=14, own cache
# sweep_rtm_p1.0.jld2. --mem=120G; the validity guard caches any overflow (|θ0|>50) as :error.

module load julia/1.12.0

export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK

julia --project=. cluster/wall_scan_cluster.jl psweep 1.0 14

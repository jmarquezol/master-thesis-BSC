#!/bin/bash
#SBATCH --job-name=wallscan_rtm_p0.3
#SBATCH --output=logs/wallscan_rtm_p0.3-%j.out
#SBATCH --error=logs/wallscan_rtm_p0.3-%j.err
#SBATCH --account=bsc21
#SBATCH --qos=gp_bsccase
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16

# p-sweep arm: p=0.3, RTM route, warm-started and checkpointed, all the way to T=20. See
# submit_rtm_p0.0.slurm for the full explanation — same design, just a different NNN coupling.
# NOTE: block_transfer_eigs_adaptive will auto-escalate k=4->6 here whenever needed (NB5's finding
# that larger p can produce a k=4 block with no tower member besides lambda0) — expect this arm to
# be more expensive per-T than p=0.0/0.1 for that reason, on top of T reaching further (20 vs 14).

cd "${SLURM_SUBMIT_DIR:-.}" || exit 1
mkdir -p logs

module load julia/1.12.0

export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK

julia --project=. cluster/wall_scan_cluster.jl psweep 0.3 20

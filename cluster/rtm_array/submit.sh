#!/bin/bash
# Submits the workers and dependent collector for the RTM wall scan.

# Move to the script's directory so it can be launched from anywhere
cd "$(dirname "$0")" || exit 1

# Create logs directory if it doesn't exist
mkdir -p logs

JOB_IDS=""

for TASK_ID in {1..13}; do
  HOURS=$((7 + TASK_ID))
  HOURS_STR=$(printf "%02d" $HOURS)
  
  JOB_ID=$(sbatch --parsable <<EOF
#!/bin/bash
#SBATCH --job-name=wallscan_rtm_${TASK_ID}
#SBATCH --output=logs/wallscan_rtm-${TASK_ID}-%j.out
#SBATCH --error=logs/wallscan_rtm-${TASK_ID}-%j.err
#SBATCH --account=bsc21
#SBATCH --qos=gp_bsccase
#SBATCH --time=${HOURS_STR}:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16

module purge
module load julia/1.12.0

export JULIA_NUM_THREADS=\$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=\$SLURM_CPUS_PER_TASK

julia --project=../.. worker.jl ${TASK_ID}
EOF
  )
  
  echo "Submitted worker task ${TASK_ID} with ${HOURS_STR}h (Job ID: $JOB_ID)"
  
  if [ -z "$JOB_IDS" ]; then
    JOB_IDS="$JOB_ID"
  else
    JOB_IDS="${JOB_IDS}:${JOB_ID}"
  fi
done

COLLECT_JOB=$(sbatch --parsable --dependency=afterany:$JOB_IDS <<EOF
#!/bin/bash
#SBATCH --job-name=wallscan_rtm_collect
#SBATCH --output=logs/wallscan_rtm_collect-%j.out
#SBATCH --error=logs/wallscan_rtm_collect-%j.err
#SBATCH --account=bsc21
#SBATCH --qos=gp_bsccase
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

module purge
module load julia/1.12.0

julia --project=../.. collect.jl
EOF
)

echo "Submitted collector job: $COLLECT_JOB (dependent on $JOB_IDS)"

#!/bin/bash
#SBATCH --job-name=hsmm_collect
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=logs/collect_%j.out
#SBATCH --error=logs/collect_%j.err

set -euo pipefail

export THESIS_DIR=$HOME/qREML
export SCENARIO_ID=${SCENARIO_ID:-1}
export SIM_RUNS=100

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate hsmm

Rscript "$THESIS_DIR/collect_runs.R"

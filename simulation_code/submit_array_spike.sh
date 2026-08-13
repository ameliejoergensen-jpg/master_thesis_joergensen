#!/bin/bash
#SBATCH --job-name=hsmm_array
#SBATCH --array=1-100
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=logs/run_%A_%a.out
#SBATCH --error=logs/run_%A_%a.err

set -euo pipefail

export THESIS_DIR=$HOME/qREML
export SCENARIO_ID=${SCENARIO_ID:-1}
export SIM_RUNS=100
export N_SELECT=100

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate hsmm

Rscript "$THESIS_DIR/git_run_array_spike.R"

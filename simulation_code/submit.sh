#!/bin/bash

set -euo pipefail

cd "$HOME/qREML"
mkdir -p logs sim_results
export SCENARIO_ID=${SCENARIO_ID:-1}

jid=$(sbatch --parsable submit_array.sh)
echo "array job: $jid"

cid=$(sbatch --parsable --dependency=afterany:"$jid" submit_collect.sh)
echo "collect job: $cid"

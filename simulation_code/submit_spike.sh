#!/bin/bash
# Submit the array, then the collection job once the array has finished.
# afterany (not afterok) so the rda is still built if a few tasks fail.

set -euo pipefail

cd "$HOME/qREML"
mkdir -p logs sim_results
export SCENARIO_ID=${SCENARIO_ID:-1}

jid=$(sbatch --parsable submit_array_spike.sh)
echo "array job: $jid"

cid=$(sbatch --parsable --dependency=afterany:"$jid" submit_collect.sh)
echo "collect job: $cid"

#!/usr/bin/env bash
# Regenerate every number and figure of paper/v2.
#   bash scripts/paper/run_all.sh            # quick mode (CPU, a few minutes)
#   bash scripts/paper/run_all.sh full cuda  # sizes reported in the paper (GPU recommended)
set -euo pipefail
MODE=${1:-quick}
DEVICE=${2:-$(python -c "import torch; print('cuda' if torch.cuda.is_available() else 'cpu')")}
cd "$(dirname "$0")/../.."
if [ "$MODE" = full ] && [ ! -f results/experiments/summary.csv ]; then
  python scripts/experiment.py --dims 2,16 --channels 1,2,16,32 --device "$DEVICE" --output-dir results/experiments
fi
for s in validate_known_results inequality_search bounds_search chains numerics baselines benchmark_library; do
  echo "=== $s ($MODE, $DEVICE)"
  python "scripts/paper/$s.py" --mode "$MODE" --device "$DEVICE"
done
for ex in examples/paper/ex*.py; do echo "=== $ex"; python "$ex"; done > paper/v2/generated/examples_output.txt
echo "done: paper/v2/generated, paper/v2/figures (compile paper/v2/main.tex)"

#!/usr/bin/env bash
# Regenerate every number and figure of paper/v2.
#   bash scripts/paper/run_all.sh                     # quick mode, CPU or GPU (a few minutes)
#   bash scripts/paper/run_all.sh check cuda          # the sizes of full, minimal repetitions: tests a full run
#                                                     # (memory, all code paths) in minutes; output in .check_output/
#   bash scripts/paper/run_all.sh full cuda           # sizes reported in the paper (GPU recommended)
#   bash scripts/paper/run_all.sh full cuda tf32      # same, with TF32 tensor-core matmuls
#   ONLY="chains baselines" bash scripts/paper/run_all.sh full cuda   # a subset of the scripts
# Every script writes its own log to paper/v2/generated/logs/; a failing script does not stop the others,
# and the run ends with a summary (exit code 1 if anything failed).
set -uo pipefail
MODE=${1:-quick}
DEVICE=${2:-$(python -c "import torch; print('cuda' if torch.cuda.is_available() else 'cpu')")}
TF32=""
if [ "${3:-}" = tf32 ]; then TF32="--tf32"; fi
EXAMPLE_ENV=""
if [ "$DEVICE" = cpu ]; then EXAMPLE_ENV="CUDA_VISIBLE_DEVICES="; fi   # the examples pick CUDA when visible
# short scripts first; inequality_search and bounds_search take hours in full mode and resume after an
# interruption (finished parts are kept in generated/partial/ until the script completes)
SCRIPTS=${ONLY:-"validate_known_results numerics baselines benchmark_library chains inequality_search bounds_search"}
cd "$(dirname "$0")/../.."
OUTDIR=paper/v2/generated
OUTARGS=""
if [ "$MODE" = check ]; then OUTDIR=.check_output/generated; fi
if [ -n "$TF32" ]; then                      # a TF32 run must not overwrite the FP32 tables of the paper
  if [ "$MODE" = check ]; then OUTDIR=.check_output/generated_tf32; else OUTDIR=paper/v2/generated_tf32; fi
  OUTARGS="--out $OUTDIR --figdir $OUTDIR/figures"
fi
LOGS=$OUTDIR/logs
mkdir -p "$LOGS"
python - <<PY | tee "$LOGS/environment.txt"
import platform, torch
print("python", platform.python_version(), "| torch", torch.__version__, "| device: $DEVICE | mode: $MODE")
if torch.cuda.is_available():
    for i in range(torch.cuda.device_count()):
        p = torch.cuda.get_device_properties(i)
        print(f"cuda:{i} {p.name}, {p.total_memory / 2**30:.0f} GiB, compute capability {p.major}.{p.minor}")
    print("CUDA", torch.version.cuda, "| cuDNN", torch.backends.cudnn.version())
elif "$DEVICE".startswith("cuda"):
    raise SystemExit("device $DEVICE requested, but torch.cuda.is_available() is False")
PY
[ "${PIPESTATUS[0]}" -eq 0 ] || exit 1

if [ "$MODE" = full ] && [ ! -f results/experiments/summary.csv ]; then
  python scripts/experiment.py --dims 2,16 --channels 1,2,16,32 --device "$DEVICE" --output-dir results/experiments \
    > "$LOGS/experiment.log" 2>&1 || experiment_failed=1
fi

failed=()
if [ -n "${experiment_failed:-}" ]; then echo "experiment.py FAILED (see $LOGS/experiment.log)"; failed+=("experiment"); fi
for s in $SCRIPTS; do
  start=$(date +%s)
  echo "=== $s ($MODE, $DEVICE${TF32:+, tf32})  log: $LOGS/$s.log"
  if python "scripts/paper/$s.py" --mode "$MODE" --device "$DEVICE" $TF32 $OUTARGS > "$LOGS/$s.log" 2>&1; then
    echo "    ok in $(( $(date +%s) - start )) s"
  else
    echo "    FAILED after $(( $(date +%s) - start )) s; last lines:"; tail -n 5 "$LOGS/$s.log" | sed 's/^/      /'
    failed+=("$s")
  fi
done
if [ -z "${ONLY:-}" ]; then
  for ex in examples/paper/ex*.py; do echo "=== $ex"; env $EXAMPLE_ENV python "$ex" || failed+=("$ex"); done > "$OUTDIR/examples_output.txt" 2>&1
fi

echo
if [ ${#failed[@]} -eq 0 ]; then
  echo "done: all scripts succeeded; outputs in $OUTDIR (compile paper/v2/main.tex)"
else
  echo "done with failures: ${failed[*]} (logs in $LOGS)"; exit 1
fi

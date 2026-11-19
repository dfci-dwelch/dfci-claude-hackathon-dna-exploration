#!/usr/bin/env bash
# Full-stream validation (resumable via per-job checkpoints). Run detached:
#   setsid nohup scripts/run_full.sh > logs/full.log 2>&1 &
cd "$(dirname "$0")/.."
for l in HOP-62 HCC2935; do
  uv run python scripts/quantify.py "$l" --mode full --workers 20 >> "logs/${l}_full.log" 2>&1 && echo "full done: $l $(date)" || echo "full FAILED: $l $(date)"
done
uv run python scripts/summarize.py
echo ALL_DONE

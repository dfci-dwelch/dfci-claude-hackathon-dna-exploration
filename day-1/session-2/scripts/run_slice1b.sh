#!/usr/bin/env bash
# Determinism check (rerun region on 4 lines) then full-stream validation. Resumable.
cd "$(dirname "$0")/.."
for l in "SW 1573" HCC2935 A427 "U-2 OS"; do
  uv run python scripts/quantify.py "$l" --mode region --workers 16 > "logs/${l// /_}_region2.log" 2>&1 && echo "region done: $l" || echo "region FAILED: $l"
done
for l in HOP-62 HCC2935; do
  uv run python scripts/quantify.py "$l" --mode full --workers 16 >> "logs/${l}_full.log" 2>&1 && echo "full done: $l $(date +%T)" || echo "full FAILED: $l"
done
uv run python scripts/summarize.py; echo ALL_DONE

#!/usr/bin/env bash
# Slice 1: region mode on 5 lines, then full-stream validation on HOP-62 and HCC2935.
cd "$(dirname "$0")/.."
for l in "SW 1573" HCC2935 A427 "U-2 OS"; do
  uv run python scripts/quantify.py "$l" --mode region > "logs/${l// /_}_region.log" 2>&1 && echo "region done: $l" || echo "region FAILED: $l"
done
for l in HOP-62 HCC2935; do
  uv run python scripts/quantify.py "$l" --mode full --workers 20 > "logs/${l}_full.log" 2>&1 && echo "full done: $l" || echo "full FAILED: $l"
done
uv run python scripts/summarize.py

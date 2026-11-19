#!/usr/bin/env bash
# Slice 2: region v2 (+ per-k-mer rDNA counts) on 4 lines, then A427 full stream (independent male validation of v2).
cd "$(dirname "$0")/.."
for l in "SW 1573" HCC2935 A427 "U-2 OS"; do
  uv run python scripts/quantify.py "$l" --mode region --outname region_v2 --workers 16 > "logs/${l// /_}_region_v2.log" 2>&1 && echo "region_v2 done: $l" || echo "region_v2 FAILED: $l"
done
uv run python scripts/quantify.py A427 --mode full --workers 16 >> logs/A427_full.log 2>&1 && echo "full done: A427 $(date +%T)" || echo "full FAILED: A427"
echo ALL_DONE

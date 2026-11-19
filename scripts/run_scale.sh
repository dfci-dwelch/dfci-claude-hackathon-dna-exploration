#!/usr/bin/env bash
# Scale-up: region v2 on data/meta/selected_lines.tsv. Resumable: skips lines with stats.json; checkpoints within a line.
cd "$(dirname "$0")/.."
tail -n +2 data/meta/selected_lines.tsv | while IFS=$'\t' read -r name mid rest; do
  if [ -f "results/$mid/region_v2/stats.json" ]; then continue; fi
  SECONDS=0
  if uv run python scripts/quantify.py "$mid" --mode region --outname region_v2 --workers 16 > "logs/${mid}_region_v2.log" 2>&1; then
    echo "ok $mid $name ${SECONDS}s"
  else
    echo "FAILED $mid $name (see logs/${mid}_region_v2.log)"
  fi
done
echo ALL_DONE

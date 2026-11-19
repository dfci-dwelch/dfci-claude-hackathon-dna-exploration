#!/usr/bin/env bash
# DJ (acrocentric p-arm dosage) region mode on all selected lines. Resumable: skips lines with stats.json.
cd "$(dirname "$0")/.."
tail -n +2 data/meta/selected_lines.tsv | while IFS=$'\t' read -r name mid rest; do
  [ -f "results/$mid/dj_region/stats.json" ] && continue
  SECONDS=0
  if uv run python scripts/quantify.py "$mid" --mode region --tag dj --regions data/kmers/hg38_regions_dj.tsv \
       --outname dj_region --workers 16 > "logs/${mid}_dj_region.log" 2>&1; then
    echo "ok $mid $name ${SECONDS}s retries=$(grep -c RETRY "logs/${mid}_dj_region.log")"
  else
    echo "FAILED $mid $name (see logs/${mid}_dj_region.log)"
  fi
done
echo ALL_DONE

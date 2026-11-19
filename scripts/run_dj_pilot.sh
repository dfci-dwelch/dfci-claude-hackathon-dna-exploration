#!/usr/bin/env bash
# DJ pilot: region mode (DJ regions only) on HOP-62 and HCC2935, then a DJ full stream of HCC2935 for capture.
cd "$(dirname "$0")/.."
for l in HOP-62 HCC2935; do
  uv run python scripts/quantify.py "$l" --mode region --tag dj --regions data/kmers/hg38_regions_dj.tsv --outname dj_region --workers 16 > "logs/${l}_dj_region.log" 2>&1 && echo "dj_region done: $l" || echo "dj_region FAILED: $l"
done
uv run python scripts/quantify.py HCC2935 --mode full --tag dj --outname dj_full --workers 16 > logs/HCC2935_dj_full.log 2>&1 && echo "dj_full done: HCC2935" || echo "dj_full FAILED: HCC2935"
echo ALL_DONE

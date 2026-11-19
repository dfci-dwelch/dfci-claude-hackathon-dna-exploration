#!/usr/bin/env bash
# Download all inputs (~5 GB). Idempotent: existing complete files are skipped; partial downloads resume.
# Usage: scripts/download_data.sh [--check]     (--check: only verify what is present, download nothing)
set -euo pipefail
cd "$(dirname "$0")/.."
CHECK=${1:-}
mkdir -p data/meta data/ref data/rdna data/depmap

# fetch URL DEST : download unless DEST exists with the server's Content-Length
fetch() {
  local url=$1 dest=$2 want have
  want=$(curl -sIL "$url" | tr -d '\r' | awk 'tolower($1)=="content-length:"{v=$2} END{print v}')
  have=$( [ -f "$dest" ] && stat -c %s "$dest" || echo 0)
  if [ -n "$want" ] && [ "$have" = "$want" ]; then echo "ok      $dest ($have bytes)"; return; fi
  if [ "$CHECK" = "--check" ]; then echo "MISSING $dest (have $have, want ${want:-?})"; return; fi
  echo "get     $dest (${want:-?} bytes)"
  curl -sS -L -C - -o "$dest" "$url"
}

# 1. CCLE inventory (public bucket, no auth)
B=https://depmap-omics-ccle.s3.amazonaws.com
fetch $B/metadata/inventory.csv        data/meta/inventory.csv
fetch $B/docs/data_dictionary.json     data/meta/data_dictionary.json

# 2. References. hg38 = the Broad reference the CCLE CRAMs were aligned to (GCS copy returns 403; use AWS).
fetch https://broad-references.s3.amazonaws.com/hg38/v0/Homo_sapiens_assembly38.fasta      data/ref/Homo_sapiens_assembly38.fasta
fetch https://broad-references.s3.amazonaws.com/hg38/v0/Homo_sapiens_assembly38.fasta.fai  data/ref/Homo_sapiens_assembly38.fasta.fai
T2T=https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies
fetch $T2T/analysis_set/chm13v2.0.fa.gz              data/ref/chm13v2.0.fa.gz        # bgzipped; pysam builds .fai/.gzi
fetch $T2T/annotation/chm13v2.0_censat_v2.1.bed      data/ref/chm13v2.0_censat_v2.1.bed

# 3. rDNA units from NCBI (45S: KY962518, U13369; 5S: X12811). No Content-Length from efetch -> check non-empty FASTA.
for a in KY962518 U13369 X12811; do
  f=data/rdna/$a.fa
  if [ -s "$f" ] && head -1 "$f" | grep -q "^>$a"; then echo "ok      $f"; continue; fi
  if [ "$CHECK" = "--check" ]; then echo "MISSING $f"; continue; fi
  echo "get     $f"; curl -sS "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$a&rettype=fasta&retmode=text" -o "$f"; sleep 0.4
done

# 4. DepMap 24Q4 Public from figshare (article 27993248); the portal itself is behind a Cloudflare check.
#    Writes data/depmap/manifest.txt (name, url, md5) and verifies md5 for every file.
FILES="AvanaGuideMap.csv KYGuideMap.csv HumagneGuideMap.csv CRISPRScreenMap.csv Model.csv OmicsProfiles.csv OmicsDefaultModelProfiles.csv OmicsAbsoluteCNSegmentsProfile.csv OmicsSignatures.csv OmicsSomaticMutationsMatrixDamaging.csv OmicsSomaticMutationsMatrixHotspot.csv CRISPRGeneEffect.csv"
if [ ! -s data/depmap/manifest.txt ]; then
  curl -sS https://api.figshare.com/v2/articles/27993248 | FILES="$FILES" python3 -I -c '
import json, os, sys
want = set(os.environ["FILES"].split())
for f in json.load(sys.stdin)["files"]:
    if f["name"] in want: print(f["name"], f["download_url"], f["supplied_md5"])' > data/depmap/manifest.txt
fi
while read -r name url md5; do
  f=data/depmap/$name
  if [ -f "$f" ] && echo "$md5  $f" | md5sum -c --quiet 2>/dev/null; then echo "ok      $f (md5)"; continue; fi
  if [ "$CHECK" = "--check" ]; then echo "MISSING $f (absent or md5 mismatch)"; continue; fi
  echo "get     $f"; curl -sS -L -o "$f" "$url"
  echo "$md5  $f" | md5sum -c --quiet || { echo "MD5 MISMATCH $f"; exit 1; }
done < data/depmap/manifest.txt
echo "all inputs present"

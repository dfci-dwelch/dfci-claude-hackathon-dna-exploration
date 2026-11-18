# rDNA and satellite-array copy number in CCLE whole genomes

This project estimates 45S and 5S rDNA copy number, and the size of three chrY satellite arrays (DYZ1, DYZ2, DYZ3), in DepMap CCLE whole-genome data. It uses 31-mers that occur only in their home array in T2T-CHM13. Reads are streamed from the public `s3://depmap-omics-ccle` bucket with HTTP range requests, so no BAM or CRAM is stored locally. The results are joined to DepMap 24Q4 copy-number, TP53 and CRISPR data.

Run every command from the project root. Python dependencies are managed with **uv** (`pyproject.toml`; Python 3.13, with pysam, numba, pandas, scipy and statsmodels). `uv run` creates the environment on first use. No samtools is needed, because pysam's bundled htslib reads the remote files directly.

## Pipeline order

| Step | Command | Runtime | Output |
|---|---|---|---|
| 0. Inputs | `scripts/download_data.sh` | ~10 min (about 5 GB) | `data/ref`, `data/rdna`, `data/meta`, `data/depmap` |
| 1. Build k-mer sets | `uv run python scripts/build_kmers.py` | ~5 min, 8 cores | `data/kmers/kmers_v1.npz`, `summary_v1.tsv`, `out_hits_v1.tsv` |
| 2. Locate k-mers in hg38 | `uv run python scripts/scan_hg38.py` | ~5 min | `data/kmers/hg38_hits_v1.tsv` |
| 3. Region list v1 | `uv run python scripts/make_regions_v1.py` | seconds | `data/kmers/hg38_regions_v1.tsv` |
| 4. Control windows | `uv run python scripts/make_controls.py` | ~1 min | `data/kmers/control_windows.tsv` |
| 5. Full-stream validation lines | `uv run python scripts/quantify.py <line> --mode full` | ~8 min per 13 GB CRAM | `results/<ACH>/full/` |
| 6. Region list v2 | `uv run python scripts/make_regions_v2.py` | seconds | `data/kmers/hg38_regions_v2.tsv` |
| 7. Select lines | `uv run python scripts/select_lines.py` | ~1 min | `data/meta/selected_lines.tsv` |
| 8. Quantify all lines | `scripts/run_scale.sh` | ~3 min per line | `results/<ACH>/region_v2/` |
| 9. Summarize | `uv run python scripts/summarize.py region_v2` | seconds | `results/summary.tsv` |
| 10. DepMap covariates | `uv run python scripts/covariates.py` | ~1 min | `data/meta/covariates.tsv` |
| 11. Analysis | `uv run python scripts/analyze.py` | ~1 min | `results/analysis/*.tsv` |

Some steps depend on others. Step 6 needs full-stream results (step 5) for at least one male line; HCC2935 was used, and A427 was held out for validation. Step 11 needs the outputs of steps 9 and 10.

To redo only the analysis on existing results, run steps 9, 10 and 11.

## Scripts

### Library and diagnostics
- **`kmerlib.py`** (library). Canonical 2-bit 31-mer encoding in numpy, used to build the k-mer sets. It also contains `scan_hits`, the numba kernel that counts hits in reads: a rolling canonical k-mer with a packed-bitmap prefilter and binary search. `scan_hits` gives exactly the same hits as the numpy path and is about 100× faster.
- **`peek_header.py <s3-or-https url>`**. Prints the SAM header of a remote CRAM or BAM using a single range request. I used it to confirm that all CCLE WGS files use Broad `Homo_sapiens_assembly38`.

### K-mer sets and target regions
- **`build_kmers.py`**. Builds the five k-mer groups:
  - 45S rDNA, from KY962518;
  - 5S rDNA, from X12811;
  - DYZ1 (HSat3), DYZ2 (HSat1B) and DYZ3 (α-satellite), from the CenSat intervals on chrY.

  It then scans all of CHM13 and keeps only k-mers whose every occurrence lies in the group's home array. For 5S, the home array is chr1:227.74–228.03 Mb; CenSat doesn't annotate it, so I set the bounds from where X12811 k-mers are densest.
- **`scan_hg38.py`**. Finds where the kept k-mers occur in the Broad hg38 reference, which shows where reads carrying them are likely to map in the CRAMs.
- **`make_regions_v1.py`**. Takes the 100 kb bins holding 99.9% of each group's hg38 hits, merges adjacent bins and pads each region by one bin. New runs use v2; v1 is only needed to reproduce the v1-vs-v2 comparison.
- **`make_controls.py`**. Picks 200 × 50 kb single-copy control windows on chr1–22. Each window has no N bases within 3 Mb and stays clear of the target regions. The median read depth across these windows is the normalization baseline.
- **`make_regions_v2.py`**. Starts from region list v1 and adds the 1 Mb bins where the full streams found hits that v1 missed. Bins are added until capture reaches 99% per group. Y groups are fitted on the male full stream only, so female background doesn't drive the choice. DYZ2 needed most of the additions.

### Quantification
**`quantify.py <cell line name or ACH id>`** is the main quantifier. It has two modes:
- `--mode region` fetches the regions listed in `--regions` (default v2), the unplaced unmapped reads (`*`), and the control windows.
- `--mode full` streams every read in the file.

Other options:
- `--outname` sets the name of the results subdirectory (default: the mode name).
- `--workers` sets the number of parallel workers (default 20; I used 16).

Read rules are the same in both modes:
- Secondary, supplementary and QC-fail reads are dropped.
- **Duplicates are kept**, because Picard marks many distinct rDNA and satellite molecules as duplicates.
- Each read is counted in the chunk that contains its start position.

**The run is resumable.** Each job's result is saved to `results/<ACH>/<outname>/parts/`, and a rerun skips finished jobs. Rerun the same command after any interruption.

Outputs:
- `hits.tsv`: hits by group, contig, 1 Mb bin and mapped/unmapped.
- `coverage.tsv`: read depth in each control window.
- `kmer_counts.tsv`: per-k-mer hit counts for the rDNA groups.
- `stats.json`: run statistics.

### Summaries and analysis
- **`summarize.py [modes]`**. Converts hits to:
  - `rel_CHM13`: array size relative to CHM13v2, haploid; for chrY, relative to HG002's Y.
  - `units_per_genome` and `units_median`: rDNA units, from the mean and the median of per-k-mer counts. The median is the default for rDNA.

  Both assume that the median control depth corresponds to 2 copies; `analyze.py` corrects for ploidy. Pass a comma-separated list of modes to filter the output, e.g. `full,region_v2`.
- **`compare_modes.py <ACH> [region_dir] [regions_file]`**. Compares region-mode and full-stream hits for one line, and shows where region mode misses hits.
- **`select_lines.py`**. Picks up to 7 lines per lineage across 12 lineages, sex-balanced where possible, using a fixed random seed (7). Lines must have 24Q4 CRISPR data and a recorded sex; engineered derivatives are excluded. It always includes the 5 slice-1 lines.
- **`covariates.py`**. Builds per-line covariates from the DepMap 24Q4 absolute copy-number segments, preferring WGS profiles:
  - `P_ctrl`: the median copy number at the control windows;
  - `acro_CN`: copy number of the proximal 10 Mb of 13q, 14q, 15q, 21q and 22q, summed;
  - `cn_1q42`: copy number at the 5S array.

  It also adds TP53 status (damaging or hotspot mutation) and ploidy, WGD and aneuploidy from OmicsSignatures.
- **`analyze.py`**. Computes ploidy- and copy-number-corrected metrics:
  - `r45_acro`: 45S units per acrocentric chromosome copy;
  - `r5_1q`: 5S units per 1q42 copy;
  - `log2_5S_45S`;
  - Y arrays adjusted by `P_ctrl`.

  It then runs the sex and loss-of-Y control, the 5S–45S correlation, TP53 tests (Mann-Whitney, and OLS adjusted for lineage and sex), lineage comparisons, and a CRISPR gene-effect scan (Spearman correlation with BH-adjusted q-values). The scan is exploratory.

### Batch wrappers
- **`run_scale.sh`**. Runs region v2 on every line in `selected_lines.tsv`, skipping lines that already have `stats.json`. It writes one line per cell line to `logs/scale.log`. Start it with `scripts/run_scale.sh >> logs/scale.log 2>&1`, and check progress with `cat logs/scale.log`.
- **`run_slice1.sh`, `run_slice1b.sh`, `run_slice2.sh`, `run_full.sh`**. Historical wrappers for the validation runs: region v1 on 5 lines, the determinism reruns, the region-v2 reruns, and the full streams of HOP-62, HCC2935 and A427. Keep them for reproducibility; they aren't needed for new lines.

## Data downloads

`scripts/download_data.sh` downloads everything below into the right directories. It skips files that are already complete: it compares file size with the server's for most files, and checks the md5 sum for DepMap files. Partial downloads resume where they stopped. `scripts/download_data.sh --check` only verifies what is present and downloads nothing.

The sources, for reference:

- **`data/meta`**
  - `https://depmap-omics-ccle.s3.amazonaws.com/metadata/inventory.csv` (1 MB; one row per file; WGS rows have `datatype == "wgs"`)
- **`data/ref`**
  - `https://broad-references.s3.amazonaws.com/hg38/v0/Homo_sapiens_assembly38.fasta` and `.fai` (3.25 GB). This is the reference the CRAMs were aligned to. The Google Cloud copy returns 403.
  - `https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz` (981 MB; bgzipped; pysam creates the `.fai` and `.gzi` indexes)
  - `.../T2T/CHM13/assemblies/annotation/chm13v2.0_censat_v2.1.bed`
- **`data/rdna`**: KY962518, U13369 and X12811 from NCBI efetch, e.g. `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=KY962518&rettype=fasta`
- **`data/depmap`**: DepMap 24Q4 Public from figshare, article 27993248. The download URLs and md5 sums are in `data/depmap/manifest.txt`. The DepMap portal itself puts automated downloads behind a Cloudflare check. The figshare files are:
  - `Model.csv`
  - `OmicsProfiles.csv`
  - `OmicsDefaultModelProfiles.csv`
  - `OmicsAbsoluteCNSegmentsProfile.csv`
  - `OmicsSignatures.csv`
  - `OmicsSomaticMutationsMatrixDamaging.csv`
  - `OmicsSomaticMutationsMatrixHotspot.csv`
  - `CRISPRGeneEffect.csv`

## Validation and caveats

**Validation results.**
- Region v2 captures 99–99.9% of full-stream hits for every group on HCC2935 and on A427, the held-out line.
- Rerunning region mode on 4 lines gave identical hits and coverage.
- Female lines have Y arrays at 0.001× or less, against about 1× for males with an intact Y.

**Caveats.**
- These are cell lines with no matched normals, so results show differences between lines, not somatic changes.
- "Unique in CHM13" doesn't hold across individuals for HSat arrays; female lines show small Y-group background.
- Absolute rDNA units depend on the ploidy and copy-number correction.
- The DepMap tables (24Q4) are older than the WGS inventory.

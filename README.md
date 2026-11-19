# rDNA and satellite-array copy number in CCLE whole genomes

This project estimates 45S and 5S rDNA copy number, and the size of three chrY satellite arrays (DYZ1, DYZ2, DYZ3), in DepMap CCLE whole-genome data. It uses 31-mers that occur only in their home array in T2T-CHM13. Reads are streamed from the public `s3://depmap-omics-ccle` bucket with HTTP range requests, so no BAM or CRAM is stored locally. The results are joined to DepMap 24Q4 copy-number, TP53 and CRISPR data.

## Setup

**Requirements.** Linux or macOS, `bash`, `curl`, internet access, about 10 GB of free disk for the reference data, and ideally 8–20 cores. No AWS account or credentials are needed: the CCLE bucket is public and is read over plain HTTPS. No samtools is needed either, because pysam's bundled htslib reads the remote files directly.

**Install uv.** Python dependencies are managed with [uv](https://docs.astral.sh/uv/). Follow the [uv installation guide](https://docs.astral.sh/uv/getting-started/installation/), or use the standalone installer:

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh    # Linux / macOS
# or: pipx install uv   |   brew install uv
```

**Create the environment.** The dependencies are pinned in `pyproject.toml` and `uv.lock` (Python ≥3.13, with pysam, numba, numpy, pandas, pyarrow, scipy and statsmodels). uv downloads a suitable Python if you don't have one.

```bash
uv sync                      # optional: `uv run` also does this on first use
scripts/download_data.sh     # step 0, ~5 GB
```

Run every command from the project root.

## Quick start

- **Look at the finished results.** The [Results summary](#results-summary) gives the main findings. Everything is already in `results/` (see [Where to find things](#where-to-find-things)). Start with `results/analysis/analysis_output.txt` and `results/analysis/per_line.tsv`.
- **Redo the analysis only.** The per-line quantification results are committed, so nothing needs to be streamed from S3. On a fresh clone `data/` is empty, though, and the analysis steps read the k-mer sets, region lists, control windows and line list from it. Run steps 0–4, 6 and 7 (about 15 minutes after the download; step 6 uses the committed full-stream results), plus `build_dj_kmers.py` if you want the DJ summary. Then run steps 9 and 10, `dj_summary.py dj_region`, and step 11.
- **Quantify a new line.** Run steps 0–4 and 6 (or keep the committed region lists in `data/kmers/` if you have them), then `uv run python scripts/quantify.py "<cell line name or ACH id>" --outname region_v2`, and rerun steps 9–11. The line must be in `data/meta/inventory.csv` with `datatype == "wgs"`.

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
| 12. DJ dosage (`scripts/run_dj.sh` for all lines) | `build_dj_kmers.py`, `scan_hg38` with tag `dj`, `make_regions_v1.py dj`, then `quantify.py <line> --tag dj --regions data/kmers/hg38_regions_dj.tsv --outname dj_region`, then `dj_summary.py` | ~1 min per line | `results/dj_summary_dj_region.tsv` |

Some steps depend on others. Step 6 needs full-stream results (step 5) for at least one male line; HCC2935 was used, and A427 was held out for validation. Step 11 needs the outputs of steps 9 and 10.

To redo only the analysis on existing results, run steps 9, 10 and 11. If `data/` is missing (for example on a fresh clone), rebuild it first as described in [Quick start](#quick-start).

## Where to find things

### What is in git and what isn't
- **In git:** `scripts/`, `results/` (except the per-job `parts/` checkpoints), `logs/`, `pyproject.toml` and `uv.lock`.
- **Not in git:** all of `data/` (download it with `scripts/download_data.sh` and regenerate the rest with steps 1–7), `.venv/`, the `results/*/*/parts/` checkpoints, and the `*.zip` archives.
- Note that `data/` is gitignored as a whole, so the derived k-mer sets, region lists, control windows, selected lines and covariates also have to be rebuilt on a fresh clone.

### Final results (`results/`)
| File | Made by | Contents |
|---|---|---|
| `results/analysis/analysis_output.txt` | `analyze.py` (its printed report, saved with `uv run python scripts/analyze.py > results/analysis/analysis_output.txt`) | Human-readable summary of every test: Y-array sex check and loss-of-Y calls, rDNA medians and IQRs, DJ vs proxy, batch effects, TP53 tests, top CRISPR hits. **Read this first.** |
| `results/analysis/per_line.tsv` | `analyze.py` | One row per line (89): metadata, raw and corrected rDNA and Y metrics (`r45_dj`, `r45_acro`, `r5_1q`, `log2_5S_45S`, `DJ_abs`, `Y_*_adj`, …) and all covariates. **The main per-line table.** |
| `results/analysis/tp53_tests.tsv` | `analyze.py` | Mann-Whitney and OLS results for each metric × TP53 status (`TP53_mut`, `p53_loss`). |
| `results/analysis/by_lineage.tsv` | `analyze.py` | Median of each metric per OncoTree lineage, with line counts. |
| `results/analysis/crispr_scan.tsv` | `analyze.py` | Spearman ρ, p and BH q for every metric × CRISPR gene. Exploratory. |
| `results/summary.tsv` | `summarize.py` | One row per line × mode × group: read length, control depth, hits, `rel_CHM13` and the rDNA unit estimates (`units_best` etc.). |
| `results/dj_summary_dj_region.tsv` | `dj_summary.py` | DJ copies per line (`DJ_raw`, `DJ_gc`, `DJ_abs`) next to the proxy `acro_CN`. |
| `results/dj_summary_dj_full.tsv` | `dj_summary.py dj_full` | The same for the HCC2935 DJ full stream (validation). |

### Per-line results (`results/<ACH id>/<outname>/`)
Each directory holds `hits.tsv`, `coverage.tsv`, `kmer_counts.tsv` (newer runs only) and `stats.json`, described under `quantify.py` below. Lines are named by DepMap model ID; `data/meta/selected_lines.tsv` maps IDs to cell-line names. The `outname` directories are:

| Directory | What it is | Lines |
|---|---|---|
| `region_v2` | Region mode with region list v2. **The production run.** | all 89 |
| `dj_region` | Region mode with the DJ k-mers and region list. | all 89 |
| `full` | Full stream of every read (validation). | HOP-62, HCC2935, A427 |
| `dj_full` | DJ full stream (validation). | HCC2935 |
| `region` | Region mode with region list v1 (historical). | the 5 slice-1 lines |
| `region_run1` | First v1 run, kept for the determinism check. | 4 of the slice-1 lines |

### Intermediate data (`data/`)
| Path | Contents |
|---|---|
| `data/ref/` | Broad hg38 (`Homo_sapiens_assembly38.fasta`), CHM13v2 (`chm13v2.0.fa.gz`) and the CenSat annotation. |
| `data/rdna/` | rDNA reference units: KY962518 (45S), U13369 (older 45S reference, downloaded for comparison), X12811 (5S). |
| `data/meta/inventory.csv` | The CCLE omics file inventory: one row per file, with `cram_bam_url` and `crai_bai_url`. `quantify.py` looks lines up here. |
| `data/meta/data_dictionary.json` | Column descriptions for `inventory.csv`. |
| `data/meta/selected_lines.tsv` | The 89 analysed lines with sex, lineage and file size (step 7). |
| `data/meta/covariates.tsv` | Per-line DepMap covariates (step 10). |
| `data/kmers/kmers_<tag>.npz` | K-mer sets (`v1` = rDNA and Y groups, `dj` = DJ): k-mer, group, home and off-target counts, position in the unit. |
| `data/kmers/summary_<tag>.tsv`, `out_hits_<tag>.tsv` | K-mers kept and dropped per group, and where the dropped ones occur elsewhere in CHM13. |
| `data/kmers/hg38_hits_<tag>.tsv` | Where the kept k-mers occur in hg38 (step 2). |
| `data/kmers/hg38_regions_{v1,v2,dj}.tsv` | The fetch lists for region mode. |
| `data/kmers/control_windows.tsv` | The 200 control windows. |
| `data/depmap/` | DepMap 24Q4 tables and the CRISPR guide maps, plus `manifest.txt` with URLs and md5 sums. |
| `data/idx_cache/` | Cached `.bai`/`.crai` indexes for the remote files. Safe to delete; they are fetched again when needed. |

### Archives
The `*.zip` files in the project root are local snapshots, not tracked in git. They include the `parts/` checkpoints, so unzipping one lets `quantify.py` skip all finished jobs.
- `final-results.zip` (450 MB): the whole `results/` directory at the end of the project, including the DJ runs and `results/analysis/`. Unzip it in the project root.
- `full-results.zip` (71 MB): an earlier snapshot (2026-11-18) of the per-line directories only, before the DJ runs. Its paths have no `results/` prefix, so unzip it inside `results/`.

### Logs (`logs/`)
- `<ACH id>_region_v2.log`, `<ACH id>_dj_region.log`: per-line output of `quantify.py` from the batch wrappers, including any `RETRY` lines.
- `scale.log`, `dj_scale.log`: one status line per cell line from `run_scale.sh` and `run_dj.sh`.
- `<line>_full.log`, `<line>_region*.log`, `slice1*.log`, `slice2.log`, `full.log`, `dj_pilot.log`: the validation runs.
- `build_kmers_v1.log`, `scan_hg38_v1.log`, `summarize.log`, `dj_summary.log`: output of the corresponding steps.

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

**Index cache.** Each BAM/CRAM index is downloaded once into `data/idx_cache/`, and its local path is passed to pysam. Passing a remote index URL makes htslib save a copy into the current directory, which is how the `CDS-*.bai` files ended up in the repo; they're now gitignored.

**The run is resumable.** Each job's result is saved to `results/<ACH>/<outname>/parts/`, and a rerun skips finished jobs. Rerun the same command after any interruption. The checkpoint key covers the job's regions plus every input that changes its result: the k-mer set, the control windows, the region file, `--tag`, and `CODE_VERSION` in `quantify.py`. Changing any of them reruns the affected jobs; it never reuses stale results. Bump `CODE_VERSION` whenever the counting logic changes.

**Retries are logged.** A failed remote read (for example htslib's `truncated file`) is retried from the start of that region or window, up to 4 times. Each retry writes a `RETRY ...` line to the line's log, so `grep RETRY logs/*.log` lists them. I tested this by forcing failures mid-stream: the results matched a clean run exactly.

Outputs:
- `hits.tsv`: hits by group, contig, 1 Mb bin and mapped/unmapped.
- `coverage.tsv`: reads, aligned bases and query bases in each control window. Query bases give the read length and were added in `CODE_VERSION` 2026-11-19.3.
- `kmer_counts.tsv`: per-k-mer hit counts for the rDNA groups.
- `stats.json`: run statistics.

- **`build_dj_kmers.py`**. Builds k-mers for the distal junction (DJ), a roughly 290 kb non-satellite block that sits next to every rDNA array in CHM13 on all five acrocentrics. Kept k-mers occur only in the 10 DJ intervals. The 140,348 k-mers that occur exactly once per DJ count rDNA-bearing p-arms directly, so they don't need the q-arm copy-number proxy. Output: `data/kmers/kmers_dj.npz`. Locate the k-mers in hg38 with `uv run python -c "import sys; sys.path.insert(0,'scripts'); import scan_hg38; scan_hg38.main('dj')"`, then build the fetch list with `make_regions_v1.py dj` (9 regions, 4.4 Mb).

### Summaries and analysis
- **`summarize.py [modes]`**. Converts hits to:
  - `rel_CHM13`: array size relative to CHM13v2, haploid; for chrY, relative to HG002's Y.
  - `units_per_genome`, `units_median`, `units_gc`: rDNA units from the mean, the median, and the median over GC-matched k-mers.
  - **`units_best`, the default rDNA metric:** `units_gc` for 45S and `units_median` for 5S.

  K-mer coverage is computed per line as c_k = D/2 × (R−K+1)/R, where:
  - **R** is the line's read length. Newer runs take it from the control-window query bases. Older runs snap it to 101 bp (BAMs) or 151 bp (CRAMs), which I confirmed directly on control reads from 8 lines. Bases per read over all fetched reads is not usable, because short unplaced-unmapped reads pull it down by up to 7%. Assuming 151 bp everywhere had made the 101 bp lines about 12% too low.
  - **D** is the median depth over the control windows, after dropping the 4 windows that overlap v2 target regions.

  **Why 45S needs a GC correction.** 45S k-mer counts fall steeply with k-mer GC. On a 151 bp line, k-mers with 25 or more GC bases (of 31) get about 0.5× the hits of moderate-GC k-mers. The drop is stronger in the 101 bp BAM batch. `units_gc` therefore uses only k-mers that occur once per unit and have 11–14 GC bases of 31, matching the control windows' median GC of 38.7%. 5S shows no GC dependence, so it uses the plain median.

  All values assume the median control depth corresponds to 2 copies; `analyze.py` corrects for ploidy. Pass a comma-separated list of modes to filter the output, e.g. `full,region_v2`.
- **`dj_summary.py [outname]`**. Gives DJ copies per genome: the median per-k-mer count over the once-per-DJ k-mers that match the controls' GC, normalized like `summarize.py`. `DJ_abs` = DJ copies × `P_ctrl`/2, which should be about 10 for a normal diploid genome. It also compares the result with DepMap `acro_CN`. rDNA per DJ copy is `r45_raw / DJ_gc`; the ploidy and control-depth assumptions cancel in that ratio.
- **`compare_modes.py <ACH> [region_dir] [regions_file]`**. Compares region-mode and full-stream hits for one line. The first table is the exact read-level capture. The list of missed bins is approximate: any overlap with a region counts the whole 1 Mb bin as covered.
- **`select_lines.py`**. Picks up to 7 lines per lineage across 12 lineages, sex-balanced where possible, using a fixed random seed (7). Lines must have 24Q4 CRISPR data and a recorded sex; engineered derivatives are excluded. It always includes the 5 slice-1 lines.
- **`covariates.py`**. Builds per-line covariates from the DepMap 24Q4 absolute copy-number segments, preferring WGS profiles:
  - `P_ctrl`: the median copy number at the control windows;
  - `acro_CN`: copy number of the proximal 10 Mb of 13q, 14q, 15q, 21q and 22q, summed;
  - `cn_1q42`: copy number at the 5S array.

  It also adds:
  - ploidy, WGD and aneuploidy from OmicsSignatures;
  - `TP53_mut`, a damaging or hotspot mutation;
  - `TP53_cn` and `TP53_LOH`, copy number and LOH at the TP53 locus;
  - `p53_score` = CRISPR TP53 effect − MDM2 effect (p53 activity);
  - `p53_functional`: `p53_score` above the 95th percentile of all 765 TP53-mutant DepMap lines (threshold 0.95);
  - `p53_loss`: a mutation, a homozygous deletion, or a non-functional `p53_score`.
- **`analyze.py`**. Computes ploidy- and copy-number-corrected metrics:
  - `r45_acro`: 45S units per acrocentric chromosome copy;
  - `r5_1q`: 5S units per 1q42 copy;
  - `log2_5S_45S`;
  - Y arrays adjusted by `P_ctrl`.

  It then runs:
  - the sex and loss-of-Y control;
  - the 5S–45S correlation;
  - TP53 tests for both `TP53_mut` and `p53_loss`: Mann-Whitney, and OLS adjusted for lineage, sex and sequencing batch (101 bp vs 151 bp). It also correlates the rDNA metrics with the continuous `p53_score`;
  - the residual batch effect on each metric;
  - lineage comparisons;
  - a CRISPR gene-effect scan: Spearman correlation using each gene's own number of lines (at least 30), with BH-adjusted q-values. The scan is exploratory.

### Review checks (`scripts/review/`)
- **`s100a8_check.py`**. Tests the S100A8–5S association three ways: S100 guides against the 5S and 45S units, the association within each screen library (Avana and KY), and copy number at 1q21 and 1q42 as a confounder. It needs the guide maps from `download_data.sh`.
- **`s100a8_baseline.py`**. Tests whether S100A8 just tracks a per-line Chronos baseline, using olfactory-receptor genes, the per-line median and the per-line SD.

### Batch wrappers
- **`run_scale.sh`**. Runs region v2 on every line in `selected_lines.tsv`, skipping lines that already have `stats.json`. It writes one line per cell line to `logs/scale.log`. Start it with `scripts/run_scale.sh >> logs/scale.log 2>&1`, and check progress with `cat logs/scale.log`.
- **`run_dj.sh`**. Runs DJ region mode on every line in `selected_lines.tsv`, skipping lines that already have `results/<ACH>/dj_region/stats.json`. It writes one status line per cell line, with its retry count; start it with `scripts/run_dj.sh >> logs/dj_scale.log 2>&1`. Run `dj_summary.py` afterwards.
- **`run_dj_pilot.sh`**. DJ pilot: region mode on HOP-62 and HCC2935, plus a DJ full stream of HCC2935.
- **`run_slice1.sh`, `run_slice1b.sh`, `run_slice2.sh`, `run_full.sh`**. Historical wrappers for the validation runs: region v1 on 5 lines, the determinism reruns, the region-v2 reruns, and the full streams of HOP-62, HCC2935 and A427. Keep them for reproducibility; they aren't needed for new lines.

## Data downloads

`scripts/download_data.sh` downloads everything below into the right directories. It skips files that are already complete: it compares file size with the server's for most files, and checks the md5 sum for DepMap files. Partial downloads resume where they stopped. `scripts/download_data.sh --check` only verifies what is present and downloads nothing.

The sources, for reference:

- **`data/meta`**
  - `https://depmap-omics-ccle.s3.amazonaws.com/metadata/inventory.csv` (1 MB; one row per file; WGS rows have `datatype == "wgs"`)
  - `https://depmap-omics-ccle.s3.amazonaws.com/docs/data_dictionary.json` (column descriptions for the inventory)
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
  - `AvanaGuideMap.csv`, `KYGuideMap.csv`, `HumagneGuideMap.csv` and `CRISPRScreenMap.csv`: guide sequences and screen-to-library maps, used only by the review checks in `scripts/review/`

## Results summary

The figures below come from `results/analysis/analysis_output.txt`, `per_line.tsv` and the review checks. The analysis is exploratory. Apart from the BH-adjusted CRISPR scan, p-values aren't corrected for the many metrics and tests, so read anything near p = 0.05 as a lead, not a finding.

### Cohort
- 89 CCLE lines across 13 OncoTree lineages: mostly 7 per lineage, with 11 lung and 1 bone because the slice-1 validation lines were added.
- 57 female and 32 male lines.
- 37 lines are 101 bp BAMs and 52 are about 151 bp CRAMs.
- 85 lines have DepMap absolute copy number, so they have the CN-corrected metrics. 88 have CRISPR data.

### The method works
- **The Y arrays separate the sexes cleanly.** Female lines have at most 0.003× of an HG002 Y in every array (median ≤0.001×). Male lines with copy-number data have median DYZ3 0.71× and DYZ1 0.63×.
- **Loss of Y is common.** Of the 30 male lines with copy-number data, 14 have Y arrays below 0.2×, for example A427, NCI-H2023, KP-4 and SW579. This agrees with the known high rate of loss of Y in male cancer lines.
- **Region mode matches full streaming** (≥99% of hits on the held-out line) **and is deterministic.** See [Validation and caveats](#validation-and-caveats).

### rDNA copy number
| Metric | Median [IQR] | n | Meaning |
|---|---|---|---|
| `r45_raw` | 304 [250, 405] | 89 | 45S units per genome, assuming the control depth is 2 copies |
| `r45_dj` | 39 [29, 46] | 89 | 45S units per DJ copy, i.e. per rDNA-bearing p-arm **(preferred 45S metric)** |
| `r45_acro` | 30 [25, 41] | 85 | 45S units per acrocentric copy from the q-arm proxy |
| `r5_raw` | 177 [139, 221] | 89 | 5S units per genome |
| `r5_1q` | 82 [60, 100] | 85 | 5S units per 1q42 copy |
| `log2_5S_45S_dj` | 1.13 [0.50, 1.55] | 85 | log2 of the 5S:45S dosage ratio |

- **45S per array varies about 7-fold between lines.** `r45_dj` runs from 11 in JHOC-5 and 13 in OPM-2 to 80 in NIH:OVCAR-3 and PE/CA-PJ41; the 5th to 95th percentile spread is about 3-fold. The median `r45_raw` of about 300 units per diploid-equivalent genome is within the few-hundred range usually reported for human genomes. Absolute values are biased low by k-mer polymorphism and sequencing error (see caveats), so comparisons between lines are more reliable than the absolute numbers.
- **5S per array varies about 6-fold,** from 32 in Hs 766T to 195 in JHOC-5. JHOC-5 has both the lowest 45S and the highest 5S.
- **5S and 45S dosage are independent.** Spearman ρ is −0.02 on raw values, +0.02 after copy-number correction, and −0.04 against `r45_dj`. Nothing in these lines co-regulates the two arrays' copy number.
- **Copy-number correction matters but keeps the ranking.** Raw versus corrected ρ is 0.89 for 45S and 0.86 for 5S.

### DJ: direct p-arm count versus the q-arm proxy
- DJ copies per genome (`DJ_abs`) have a median of 10.8, against 13.0 predicted by the proxy (`acro_CN`). The two agree only moderately (ρ 0.73).
- 22 of 85 lines have fewer than 0.75× the DJ copies the proxy predicts, which indicates loss or rearrangement of the acrocentric p-arms. The largest deficits are U-2 OS and HOP-62.
- `r45_dj` correlates 0.74 with `r45_acro`, so the main conclusions hold under either normalization. `r45_dj` is preferred because it counts the arms that actually carry rDNA.

### TP53 and p53 activity
- **Binary TP53 status shows no significant rDNA difference.** This holds for every metric, for both `TP53_mut` and the broader `p53_loss`, with Mann-Whitney and with OLS adjusted for lineage, sex and batch (`tp53_tests.tsv`). The closest result is `r45_dj` with `p53_loss`: median 41.3 versus 32.7 units, MWU p = 0.075, adjusted effect +0.20 log2 (about +15%), p = 0.12.
- **The continuous p53 activity score shows a weak trend in the same direction.** Lines with more p53 activity have fewer 45S units per DJ copy (ρ = −0.22, p = 0.036, n = 88) and a higher 5S:45S ratio (ρ = +0.22, p = 0.04). The CRISPR scan agrees: TP53 gene effect against `log2_5S_45S_dj` has ρ = +0.28, p = 0.009. None of this holds for `r45_acro` or `r5_1q`, and none of it survives multiple-testing correction.
- **Reading:** at most a weak tendency for p53-deficient lines to carry more 45S per array. A study powered for this would need more lines.

### Lineage
- 45S per array differs modestly by lineage (Kruskal-Wallis p = 0.055 for `r45_dj` and p = 0.049 for `r45_acro`). CNS/Brain lines are highest (median `r45_dj` 60) and lymphoid lines lowest (28).
- 5S does not differ by lineage (p = 0.68).
- With 7 lines per lineage, and lineage confounded with sequencing batch, these are leads only. Medians per lineage are in `by_lineage.tsv`.

### Sequencing batch
After the read-length and GC corrections, the lineage-adjusted 151 bp versus 101 bp difference is not significant for `r45_dj`, `r45_acro`, `DJ_abs` or the Y arrays. A residual remains for `r5_1q`: −0.26 log2 (about −16%), p = 0.043. Batch is therefore a covariate in every adjusted model.

### CRISPR dependencies
- **Genes that would be expected to respond to rDNA dosage show nothing.** RNA Pol I subunits (POLR1A–D), UBTF, TCOF1, FBL, NCL, NOP56, MYC, Pol III subunits and GTF3C1 all have |ρ| < 0.2 and p > 0.08 against every metric. The one exception is RRN3 against `r45_acro` (ρ = −0.32, p = 0.044), but only 39 lines were screened for RRN3. Cell lines with more rDNA copies are not more or less dependent on the rRNA transcription machinery.
- **Across the genome, no gene passes q < 0.05 for any 45S metric;** the best q is about 0.07. Only **S100A8** passes for 5S: against `r5_1q`, ρ = +0.53 and q = 0.004; against `log2_5S_45S_dj`, q = 0.016. The neighbouring S100A12 and S100A4 also rank near the top for `r5_1q`.
- **The S100A8 hit survives every artifact check in `scripts/review/`:**
  - *Guide off-targeting onto rDNA:* every S100-cluster guide has at least 5 mismatches to both rDNA units.
  - *A single screen library:* the association holds within Avana (n = 75, ρ = 0.46) and within KY (n = 34, ρ = 0.63).
  - *1q copy number:* the S100 cluster sits at 1q21 and the 5S array at 1q42, but 1q21 and 1q42 copy number don't explain it. The slope barely changes after adjusting for copy number, lineage and batch.
  - *Per-line Chronos baseline:* `r5_1q` doesn't track the olfactory-receptor baseline (ρ = 0.04), and the slope stays at +0.07, p = 7e-6, after adding baseline terms.
  - *A few outlier lines:* dropping the 6 lines with the most positive S100A8 effects leaves ρ = 0.50.
- **But the effect is biologically odd.** S100A8 is not essential in these lines: its gene effects run from −0.17 to +0.09 (5th to 95th percentile), and it is mainly a myeloid gene. The association is statistically robust, but there is no mechanism for it, so treat it as an unexplained hypothesis, not a dependency.

### Bottom line
CCLE lines differ several-fold in 45S and 5S dosage per array, and the two arrays vary independently. Neither tracks TP53 status strongly or predicts dependency on the rRNA machinery. The direct DJ count shows that the q-arm copy-number proxy often misses p-arm loss. The robust but unexplained S100A8–5S association is the one genome-wide CRISPR signal.

## Validation and caveats

**Validation results.**
- Region v2 captures 99–99.9% of full-stream hits for every group on HCC2935 and on A427, the held-out line. The one exception is HCC2935 DYZ2 at 98.95%; HCC2935 is the line the v2 bins were fitted on. These figures are read-level, from `compare_modes.py`.
- Rerunning region mode on 4 lines gave identical hits and coverage.
- Female lines have Y arrays at 0.001× or less, against about 1× for males with an intact Y.

**Caveats.**
- These are cell lines with no matched normals, so results show differences between lines, not somatic changes.
- "Unique in CHM13" doesn't hold across individuals for HSat arrays; female lines show small Y-group background.
- Absolute rDNA units depend on the ploidy and copy-number correction. That correction uses a proxy, proximal q-arm copy number, which can't see p-arm losses or translocations. DJ has now been run on all 89 lines, and `r45_dj` (45S per DJ copy) is the preferred 45S metric.
  - DJ agrees with the proxy only moderately (Spearman 0.73). 22 of 85 lines have fewer than 0.75× the DJ copies the proxy predicts, which means p-arm loss.
  - The biggest deficits are U-2 OS (4.9 DJ copies vs 15.2 from the proxy) and HOP-62 (5.0 vs 13.3). For those lines, proxy-based per-arm rDNA is about 3× too low.
  - On 11 near-diploid lines DJ has a median of 9.4 copies, against the expected 10.
- **Sequencing batch is confounded with lineage.** 37 lines are 101 bp BAMs and 52 are about 151 bp CRAMs, and the mix differs by lineage. The read-length and GC corrections remove most of the difference, but a lineage-adjusted residual remains for 5S (about −16% in the 151 bp batch, p≈0.05), so batch is a covariate in the adjusted models.
- Y arrays get no GC correction, since no per-k-mer counts are saved for them.
- No correction is applied for sequencing errors or k-mer polymorphism, which lower all absolute values.
- The DepMap tables (24Q4) are older than the WGS inventory.

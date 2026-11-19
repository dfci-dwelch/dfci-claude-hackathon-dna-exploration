"""Count array-specific k-mer hits in a remote CCLE WGS CRAM/BAM.

  mode=region : fetch only hg38 target regions (--regions, default data/kmers/hg38_regions_v2.tsv), unplaced
                unmapped reads ('*'), and single-copy control windows.
  mode=full   : stream every read (genome split into chunks, plus all small contigs, plus '*').

Read rules (identical in both modes): skip secondary/supplementary/QC-fail; duplicates are KEPT
(Picard flags many distinct rDNA/satellite molecules as duplicates because copies are identical).
A read is attributed to the chunk containing its reference_start, so overlapping fetches never
double count. Placed-unmapped reads (mate mapped) are attributed to their mate's position.

Outputs results/<model_id>/<mode>/: hits.tsv (group, contig, bin_1Mb, unmapped, hits),
coverage.tsv (control window reads, aligned bases, query bases), kmer_counts.tsv (per-k-mer hits for rDNA and DJ groups), stats.json.
Output dir name = --outname (default: mode).
"""
import sys, os, json, time, argparse, multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor
import numpy as np, pandas as pd, pysam
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kmerlib import encode, scan_hits, make_packed_bitmap

HG38 = "data/ref/Homo_sapiens_assembly38.fasta"
HTTP = "https://depmap-omics-ccle.s3.amazonaws.com/"
BITS, BATCH_BASES, FULL_CHUNK = 27, 6_000_000, 20_000_000
CODE_VERSION = "2026-11-19.3"  # .2: per-k-mer counts for DJ; .3: control query_bases (read length)  # bump when counting logic changes: invalidates checkpoints
SKIP = 0x100 | 0x800 | 0x200

# globals inherited by forked workers
T = G = BM = None; NAMES = []; CONTROLS = None; PERK = None  # PERK: bool mask of targets tracked per k-mer


def load_targets(tag):
    global T, G, BM, NAMES
    d = np.load(f"data/kmers/kmers_{tag}.npz")
    keep = (d["n_out"] == 0) & (d["n_home"] > 0)
    T, G, NAMES = d["kmer"][keep], d["group"][keep], list(d["group_names"])
    BM = make_packed_bitmap(T, BITS)
    global PERK
    PERK = np.array([NAMES[g].startswith("rDNA") or NAMES[g] == "DJ" for g in G])  # groups with per-k-mer counts
    return d


IDX_CACHE = "data/idx_cache"


def local_index(row):
    """Download the BAI/CRAI once into data/idx_cache/ and return its path. Passing a local index stops htslib
    from saving a copy of the remote index into the current working directory."""
    import urllib.request
    url = row["crai_bai_url"].replace("s3://depmap-omics-ccle/", HTTP)
    dest = os.path.join(IDX_CACHE, os.path.basename(url))
    if not os.path.exists(dest):
        os.makedirs(IDX_CACHE, exist_ok=True)
        tmp = f"{dest}.{os.getpid()}.tmp"
        urllib.request.urlretrieve(url, tmp)
        os.replace(tmp, dest)
    return dest


def open_aln(row):
    u = row["cram_bam_url"].replace("s3://depmap-omics-ccle/", HTTP)
    i = row.get("_local_index") or local_index(row)
    for attempt in range(5):
        try:
            return pysam.AlignmentFile(u, index_filename=i, reference_filename=HG38)
        except OSError as e:
            err = e; log_retry("open", row["model_id"], attempt, e)
            time.sleep(5 * (attempt + 1))
    raise err


def log_retry(what, where, attempt, err):
    print(f"RETRY {what} {where} attempt={attempt + 1} error={err!r}", file=sys.stderr, flush=True)


class Acc:
    """Accumulates reads into a concatenated buffer and tallies k-mer hits per (group, contig, bin, unmapped)."""
    def __init__(self):
        self.buf, self.meta, self.nb = [], [], 0
        self.tally = {}; self.reads = 0; self.bases = 0

    def add(self, seq, contig, pos, unmapped):
        self.buf.append(seq); self.meta.append((contig, pos // 1_000_000 if pos >= 0 else -1, unmapped))
        self.nb += len(seq) + 1
        if self.nb >= BATCH_BASES: self.flush()

    def flush(self):
        if not self.buf: return
        s = "N".join(self.buf); self.reads += len(self.buf); self.bases += len(s) - len(self.buf) + 1
        codes = encode(s)
        Phit, idx = scan_hits(codes, T, BM, BITS)  # exact; identical output to the earlier numpy path
        if len(idx):
            starts = np.cumsum([0] + [len(x) + 1 for x in self.buf[:-1]])
            ridx = np.searchsorted(starts, Phit, side="right") - 1
            u, cnt = np.unique(np.stack([ridx, G[idx].astype(np.int64)]), axis=1, return_counts=True)
            for r, g, n in zip(u[0], u[1], cnt):
                key = (NAMES[g],) + self.meta[r]
                self.tally[key] = self.tally.get(key, 0) + int(n)
            pk = idx[PERK[idx]]
            if len(pk):
                ui, uc = np.unique(pk, return_counts=True)
                for i, n in zip(ui.tolist(), uc.tolist()):
                    key = ("__kmer__", i); self.tally[key] = self.tally.get(key, 0) + n
        self.buf, self.meta, self.nb = [], [], 0


def work(args):
    row, units, with_controls = args
    t0 = time.time(); f = open_aln(row); acc = Acc(); cov = []
    for contig, s, e in units:
        for attempt in range(4):
            try:
                it = f.fetch(contig="*") if contig == "*" else f.fetch(contig, s, e)
                sub = Acc()
                for r in it:
                    if r.flag & SKIP: continue
                    p = r.reference_start
                    if contig != "*" and not (s <= p < e): continue
                    q = r.query_sequence
                    if q: sub.add(q, contig if contig != "*" else "*", p if contig != "*" else -1, r.is_unmapped)
                sub.flush(); break
            except OSError as err:  # not 'e': that is the region end
                log_retry("fetch", f"{row['model_id']}:{contig}:{s}-{e}", attempt, err)
                if attempt == 3: raise
                time.sleep(10 * (attempt + 1)); f = open_aln(row)
        for k, v in sub.tally.items(): acc.tally[k] = acc.tally.get(k, 0) + v
        acc.reads += sub.reads; acc.bases += sub.bases
    if with_controls:
        for _, (c, s, e) in CONTROLS.iterrows():
            for attempt in range(4):  # transient 'truncated file' errors on remote reads: retry the window from scratch
                try:
                    b = n = qb = 0
                    for r in f.fetch(c, s, e):
                        if r.flag & SKIP or r.is_unmapped or not (s <= r.reference_start < e): continue
                        b += r.query_alignment_length; n += 1; qb += r.query_length
                    break
                except OSError as err:
                    log_retry("control", f"{row['model_id']}:{c}:{s}-{e}", attempt, err)
                    if attempt == 3: raise
                    time.sleep(10 * (attempt + 1)); f = open_aln(row)
            cov.append((c, s, e, n, b, qb))
    return acc.tally, acc.reads, acc.bases, cov, time.time() - t0


def full_units(f):
    big = [(c, L) for c, L in zip(f.references, f.lengths) if L >= 5_000_000]
    units = [[(c, s, min(s + FULL_CHUNK, L))] for c, L in big for s in range(0, L, FULL_CHUNK)]
    small = [(c, 0, L) for c, L in zip(f.references, f.lengths) if L < 5_000_000]
    for i in range(8): units.append(small[i::8])
    units.append([("*", 0, 0)])
    return units


def region_units(path):
    r = pd.read_csv(path, sep="\t").sort_values(["contig", "start"])
    merged = []
    for c, s, e in r[["contig", "start", "end"]].itertuples(index=False):
        if merged and merged[-1][0] == c and s <= merged[-1][2]: merged[-1][2] = max(merged[-1][2], e)
        else: merged.append([c, s, e])
    units = [[tuple(m)] for m in merged]
    units.append([("*", 0, 0)])
    return units


def main():
    global CONTROLS
    ap = argparse.ArgumentParser()
    ap.add_argument("line"); ap.add_argument("--mode", choices=["region", "full"], default="region")
    ap.add_argument("--workers", type=int, default=20); ap.add_argument("--tag", default="v1")
    ap.add_argument("--regions", default="data/kmers/hg38_regions_v2.tsv"); ap.add_argument("--outname", default=None)
    a = ap.parse_args()
    inv = pd.read_csv("data/meta/inventory.csv"); inv = inv[inv.datatype == "wgs"]
    row = inv[(inv.cell_line_name == a.line) | (inv.model_id == a.line)].iloc[0].to_dict()
    row["_local_index"] = local_index(row)  # fetched once here, before workers fork
    load_targets(a.tag)
    CONTROLS = pd.read_csv("data/kmers/control_windows.tsv", sep="\t")
    out = f"results/{row['model_id']}/{a.outname or a.mode}"; os.makedirs(out, exist_ok=True)
    f = open_aln(row)
    units = full_units(f) if a.mode == "full" else region_units(a.regions)
    # contig_lengths for clipping region ends
    L = dict(zip(f.references, f.lengths))
    units = [[(c, s, min(e, L.get(c, e))) for c, s, e in u] for u in units]
    jobs = [(row, u, False) for u in units]
    # controls in their own job, unless full mode covers them anyway (still measured identically)
    jobs.append((row, [], True))
    print(f"{row['cell_line_name']} {row['model_id']} {a.mode}: {len(jobs)} jobs", flush=True)
    tally, reads, bases, cov, t0 = {}, 0, 0, [], time.time()
    # per-job checkpoints: a restarted run reuses finished jobs. The key covers the job's units plus every
    # input that changes its result (k-mer set, control windows, region file, tag, code version), so editing
    # any of them can never silently reuse stale parts.
    import pickle, hashlib
    parts = f"{out}/parts"; os.makedirs(parts, exist_ok=True)
    def fhash(path): return hashlib.md5(open(path, "rb").read()).hexdigest()
    fp = (CODE_VERSION, a.tag, a.mode, fhash(f"data/kmers/kmers_{a.tag}.npz"), fhash("data/kmers/control_windows.tsv"),
          fhash(a.regions) if a.mode == "region" else None)
    def ppath(i):
        h = hashlib.md5(repr((jobs[i][1], jobs[i][2], fp)).encode()).hexdigest()[:10]
        return f"{parts}/{i:04d}_{h}.pkl"
    results = {}
    for i in range(len(jobs)):
        if os.path.exists(ppath(i)): results[i] = pickle.load(open(ppath(i), "rb"))
    todo = [i for i in range(len(jobs)) if i not in results]
    print(f"  reusing {len(results)} checkpointed jobs, {len(todo)} to run", flush=True)
    from concurrent.futures import as_completed
    with ProcessPoolExecutor(a.workers, mp_context=mp.get_context("fork")) as ex:
        futs = {ex.submit(work, jobs[i]): i for i in todo}
        for nd, fu in enumerate(as_completed(futs), 1):
            i = futs[fu]; res = fu.result(); results[i] = res
            tmp = ppath(i) + ".tmp"; pickle.dump(res, open(tmp, "wb")); os.replace(tmp, ppath(i))
            if nd % 10 == 0 or a.mode == "region" or nd == len(todo):
                done_reads = sum(r[1] for r in results.values())
                print(f"  {len(results)}/{len(jobs)} reads={done_reads:,} elapsed={time.time()-t0:.0f}s", flush=True)
    for i in range(len(jobs)):
        tl, n, b, c, dt = results[i]
        for k, v in tl.items(): tally[k] = tally.get(k, 0) + v
        reads += n; bases += b; cov += c
    kc = {k[1]: v for k, v in tally.items() if k[0] == "__kmer__"}
    tally = {k: v for k, v in tally.items() if k[0] != "__kmer__"}
    pk_idx = np.nonzero(PERK)[0]
    pd.DataFrame({"kmer_index": pk_idx, "group": [NAMES[G[i]] for i in pk_idx], "hits": [kc.get(int(i), 0) for i in pk_idx]}) \
        .to_csv(f"{out}/kmer_counts.tsv", sep="\t", index=False)
    pd.DataFrame([k + (v,) for k, v in tally.items()], columns=["group", "contig", "bin_Mb", "unmapped", "hits"]) \
        .sort_values("hits", ascending=False).to_csv(f"{out}/hits.tsv", sep="\t", index=False)
    pd.DataFrame(cov, columns=["contig", "start", "end", "reads", "aligned_bases", "query_bases"]).to_csv(f"{out}/coverage.tsv", sep="\t", index=False)
    json.dump({"line": row["cell_line_name"], "model_id": row["model_id"], "sex": row["sex"], "mode": a.mode, "regions": a.regions if a.mode == "region" else None,
               "reads": reads, "bases": bases, "seconds": round(time.time() - t0), "code_version": CODE_VERSION, "input_fingerprint": list(fp)}, open(f"{out}/stats.json", "w"), indent=1)
    print("done", out, f"{time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()

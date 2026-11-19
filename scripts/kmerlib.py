"""Vectorized canonical k-mer encoding (k<=31) with numpy.

A k-mer is packed 2 bits/base into uint64 (A=0,C=1,G=2,T=3); canonical = min(fwd, revcomp).
Windows containing any non-ACGT base are dropped.
"""
import numpy as np

K = 31
_LUT = np.full(256, 4, dtype=np.uint8)
for i, b in enumerate(b"ACGT"):
    _LUT[b] = i
    _LUT[b + 32] = i  # lowercase (soft-masked)


def encode(seq) -> np.ndarray:
    """bytes/str -> uint8 codes, 4 = invalid."""
    if isinstance(seq, str):
        seq = seq.encode()
    return _LUT[np.frombuffer(seq, dtype=np.uint8)]


def canonical_kmers(codes: np.ndarray, k: int = K):
    """Return (kmers uint64, start positions int64) for all valid windows of `codes`."""
    n = len(codes) - k + 1
    if n <= 0:
        return np.empty(0, np.uint64), np.empty(0, np.int64)
    bad = (codes == 4).astype(np.int32)
    cs = np.concatenate([[0], np.cumsum(bad)])
    valid = (cs[k:] - cs[:-k]) == 0
    c = np.where(codes == 4, 0, codes).astype(np.uint64)
    fwd = np.zeros(n, np.uint64)
    rev = np.zeros(n, np.uint64)
    three = np.uint64(3)
    for j in range(k):
        s = c[j:j + n]
        fwd = (fwd << np.uint64(2)) | s
        rev |= (three - s) << np.uint64(2 * j)
    can = np.minimum(fwd, rev)
    pos = np.nonzero(valid)[0]
    return can[pos], pos.astype(np.int64)


def lookup(sorted_targets: np.ndarray, kmers: np.ndarray):
    """Return (mask of kmers present in targets, their indices into targets)."""
    if len(sorted_targets) == 0 or len(kmers) == 0:
        return np.zeros(len(kmers), bool), np.empty(0, np.int64)
    idx = np.searchsorted(sorted_targets, kmers)
    idx[idx == len(sorted_targets)] = 0
    hit = sorted_targets[idx] == kmers
    return hit, idx[hit]


def canonical_kmers_strided(codes: np.ndarray, stride: int = 1, k: int = K):
    """Like canonical_kmers but only evaluates windows starting at positions p % stride == 0.
    Unbiased (x stride) estimator of hit totals when reads are concatenated at arbitrary offsets."""
    n = len(codes) - k + 1
    if n <= 0:
        return np.empty(0, np.uint64), np.empty(0, np.int64)
    bad = (codes == 4).astype(np.int32)
    cs = np.concatenate([[0], np.cumsum(bad)])
    P = np.arange(0, n, stride)
    P = P[(cs[P + k] - cs[P]) == 0]
    c = codes.astype(np.uint64)
    fwd = np.zeros(len(P), np.uint64); rev = np.zeros(len(P), np.uint64)
    three = np.uint64(3)
    for j in range(k):
        s = c[P + j]
        fwd = (fwd << np.uint64(2)) | s
        rev |= (three - s) << np.uint64(2 * j)
    return np.minimum(fwd, rev), P


_MULT = np.uint64(0x9E3779B97F4A7C15)


def bitmap_hash(kmers: np.ndarray, bits: int) -> np.ndarray:
    return ((kmers * _MULT) >> np.uint64(64 - bits)).astype(np.int64)


# ---- numba single-pass kernel (rolling canonical k-mer + packed-bitmap prefilter + binary search) ----
import numba as nb


@nb.njit(cache=True, nogil=True)
def scan_hits(codes, T, packed, bits, k=31):
    """Return (positions, target indices) of every window in `codes` whose canonical k-mer is in sorted T.
    Windows containing code 4 (non-ACGT / read separator) are skipped. Exact (no subsampling)."""
    n = len(codes)
    out_p = np.empty(max(n // 4, 1024), np.int64); out_i = np.empty(max(n // 4, 1024), np.int64); m = 0
    mask = (np.uint64(1) << np.uint64(2 * k)) - np.uint64(1)
    shift = np.uint64(2 * (k - 1)); mult = np.uint64(0x9E3779B97F4A7C15); hb = np.uint64(64 - bits)
    fwd = np.uint64(0); rev = np.uint64(0); run = 0
    nT = len(T)
    for i in range(n):
        c = codes[i]
        if c > 3:
            run = 0; fwd = np.uint64(0); rev = np.uint64(0); continue
        cu = np.uint64(c)
        fwd = ((fwd << np.uint64(2)) | cu) & mask
        rev = (rev >> np.uint64(2)) | ((np.uint64(3) - cu) << shift)
        run += 1
        if run < k:
            continue
        can = fwd if fwd < rev else rev
        h = (can * mult) >> hb
        if (packed[h >> np.uint64(3)] >> (h & np.uint64(7))) & np.uint8(1):
            lo = 0; hi = nT
            while lo < hi:
                mid = (lo + hi) >> 1
                if T[mid] < can: lo = mid + 1
                else: hi = mid
            if lo < nT and T[lo] == can:
                if m == len(out_p):
                    out_p = np.concatenate((out_p, np.empty(len(out_p), np.int64)))
                    out_i = np.concatenate((out_i, np.empty(len(out_i), np.int64)))
                out_p[m] = i - k + 1; out_i[m] = lo; m += 1
    return out_p[:m], out_i[:m]


def make_packed_bitmap(T, bits):
    bm = np.zeros(1 << bits, bool); bm[bitmap_hash(T, bits)] = True
    return np.packbits(bm, bitorder="little")

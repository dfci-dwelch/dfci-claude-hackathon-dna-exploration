"""Fetch the first bytes of a remote CRAM/BAM via HTTP range and print the SAM header text."""
import sys, urllib.request, zlib, gzip, io, re
url = sys.argv[1].replace("s3://depmap-omics-ccle/", "https://depmap-omics-ccle.s3.amazonaws.com/")
req = urllib.request.Request(url, headers={"Range": "bytes=0-1048575"})
buf = urllib.request.urlopen(req).read()
if buf[:4] == b"CRAM":
    print("CRAM version", buf[4], buf[5])
    # header block is usually raw or gzip inside first container; scan for gzip magic or raw '@HD'
    i = buf.find(b"@HD")
    if i >= 0:
        txt = buf[i:].split(b"\x00")[0].decode(errors="replace")
    else:
        j = buf.find(b"\x1f\x8b", 26)
        txt = zlib.decompressobj(31).decompress(buf[j:]).decode(errors="replace")
else:
    out = b""; d = io.BytesIO(buf)
    try:
        with gzip.GzipFile(fileobj=d) as g:
            out = g.read(4_000_000)
    except (EOFError, OSError):
        pass
    l = int.from_bytes(out[4:8], "little"); txt = out[8:8+l].decode(errors="replace")
lines = txt.splitlines()
sq = [l for l in lines if l.startswith("@SQ")]
print("n @SQ:", len(sq))
for l in lines[:3] + sq[:3] + [s for s in sq if "GL000220" in s or "chrEBV" in s or "HLA" in s][:3] + sq[-2:]:
    print(l[:300])
print("decoy contigs:", sum("decoy" in s for s in sq), " HLA:", sum("HLA-" in s for s in sq), " alt:", sum("_alt" in s for s in sq))
for l in lines:
    if l.startswith("@PG") and ("bwa" in l.lower() or "dragen" in l.lower()): print(l[:400]); break

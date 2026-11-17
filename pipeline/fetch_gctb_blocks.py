#!/usr/bin/env python3
"""Fetch only the GCTB FullLDM members the 250 loci need, straight out of the server's zip.

The UKB EUR 13M FullLDM is one 206.6 GB zip on gctbhub. The server answers HTTP range
requests, and every member is compressed on its own, so `zipfile` over a seekable
range-request reader extracts single members without downloading the archive: the index is
4 requests, and each block streams out on its own. zipfile checks each member's CRC as it
reads.

Which blocks: every block in `../inputs/ref_b37_1588blocks.pos` (the server's own GRCh37
block map) that overlaps a locus window. For the 250 loci that is 229 blocks, 39.8 GB
compressed and 51.1 GB unpacked. `snp.info` and `ldm.info` always come too.

Resumable: a member already on disk at its full uncompressed size is skipped.

Usage:
  python3 fetch_gctb_blocks.py --loci ../inputs/loci_250.txt --out <dir> [--list-only]
"""
import argparse
import io
import os
import urllib.request
import zipfile

import pandas as pd

URL = ("https://gctbhub.cloud.edu.au/data/SBayesRC/resources/GWFM/LD/Imputed13M/"
       "ukbEUR_13M_FullLDM.zip")
HERE = os.path.dirname(os.path.abspath(__file__))


class HTTPRange(io.RawIOBase):
    """Seekable read-only file over HTTP range requests."""

    def __init__(self, url):
        self.url, self.pos = url, 0
        head = urllib.request.urlopen(urllib.request.Request(url, method="HEAD"))
        self.size = int(head.headers["Content-Length"])

    def seekable(self):
        return True

    def readable(self):
        return True

    def tell(self):
        return self.pos

    def seek(self, off, whence=0):
        self.pos = {0: off, 1: self.pos + off, 2: self.size + off}[whence]
        return self.pos

    def readinto(self, b):
        if self.pos >= self.size or not len(b):
            return 0
        end = min(self.pos + len(b), self.size) - 1
        req = urllib.request.Request(self.url, headers={"Range": f"bytes={self.pos}-{end}"})
        data = urllib.request.urlopen(req).read()
        b[:len(data)] = data
        self.pos += len(data)
        return len(data)


def blocks_for(loci_file, pos_file):
    b = pd.read_csv(pos_file, sep=r"\s+")
    need = set()
    for line in open(loci_file):
        loc = line.strip()
        if not loc:
            continue
        c, s, e = loc.split("_")
        c, s, e = int(c.replace("chr", "")), int(s), int(e)
        need |= set(b[(b.Chrom == c) & (b.StartBP <= e) & (b.EndBP >= s)].Block.astype(int))
    return sorted(need)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--loci", default=os.path.join(HERE, "..", "inputs", "loci_250.txt"))
    ap.add_argument("--pos", default=os.path.join(HERE, "..", "inputs", "ref_b37_1588blocks.pos"))
    ap.add_argument("--out", required=True)
    ap.add_argument("--list-only", action="store_true")
    args = ap.parse_args()

    blocks = blocks_for(args.loci, args.pos)
    z = zipfile.ZipFile(io.BufferedReader(HTTPRange(URL), buffer_size=16 << 20))
    info = {i.filename: i for i in z.infolist()}
    names = ["ldm13M/ldm.info", "ldm13M/snp.info"] + [f"ldm13M/block{k}.ldm.bin" for k in blocks]
    comp = sum(info[n].compress_size for n in names)
    raw = sum(info[n].file_size for n in names)
    print(f"{len(blocks)} blocks + 2 info files: {comp / 1e9:.1f} GB compressed, "
          f"{raw / 1e9:.1f} GB unpacked")
    if args.list_only:
        print(" ".join(str(k) for k in blocks))
        return
    os.makedirs(args.out, exist_ok=True)
    for i, n in enumerate(names, 1):
        dst = os.path.join(args.out, os.path.basename(n))
        if os.path.exists(dst) and os.path.getsize(dst) == info[n].file_size:
            continue
        with z.open(n) as src, open(dst + ".part", "wb") as fh:
            while True:
                chunk = src.read(16 << 20)
                if not chunk:
                    break
                fh.write(chunk)
        os.replace(dst + ".part", dst)
        print(f"[{i}/{len(names)}] {n} {info[n].file_size / 1e6:.0f} MB")


if __name__ == "__main__":
    main()

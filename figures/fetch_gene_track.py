#!/usr/bin/env python3
"""Fetch the Fig. 2 gene tracks from Ensembl (GRCh38) into results/fig2_loci/fig2_<locus>_genes.csv.

These tables were originally pulled by hand; this makes
them reproducible. Two things the hand-pulled versions did not do:

  * ONLY NAMED GENES. Ensembl returns novel lncRNAs and processed pseudogenes with no
    `external_name`, and the track was labelling them with raw ENSG accessions. An ENSG
    number tells a reader nothing, so unnamed features are dropped.

  * THE RIGHT GENOME BUILD. The PGC3-SCZ sumstats, the locus boundaries and every bp in
    this pipeline are GRCh37; the hand-pulled tables came from Ensembl's default REST
    endpoint, which is GRCh38, and were wrong at all three loci. chr12's lead rs4766428
    sits at 110,723,245 in GRCh37 (inside ATP2A2) and 110,285,440 in GRCh38, so the old
    track put it in PPP1CC, 372 kb away. `--assembly` now selects the endpoint and the
    fetch REFUSES to write unless the locus lead's position on that assembly equals the
    bp our own table carries.

  * NO CLONE-NAME OR PSEUDOGENE CLUTTER. GRCh37's annotation is full of clone-based
    placeholders (RP11-478C19.2, AC144522.1) and processed pseudogenes (RPL31P49,
    RN7SL769P). They carry no meaning for a reader and they collided with the real labels
    -- chr12 alone had 38 features, 16 of them noise. Protein-coding genes are always
    kept; other biotypes only when the symbol is a real one.

  * FLANKING GENES. chr1:238,935,959-239,252,559 contains no named gene at all -- it is a
    desert, nearest protein-coding CHRM3 134 kb past the right edge -- which would leave an
    empty panel. Genes within --flank-kb of the interval are kept, marked `flank=1`, and the
    plotter draws them clipped to the axis edge with their distance.

Usage:
  python3 fetch_gene_track.py --locus chr1_238935959_239252559 [...]
"""
import argparse
import json
import os
import re
import urllib.request

import pandas as pd

HOST = {"grch37": "https://grch37.rest.ensembl.org",
        "grch38": "https://rest.ensembl.org"}
# clone-based placeholder symbols: RP11-478C19.2, AC144522.1, CTD-2008L17.1, AP001267.3
CLONE = re.compile(r"^(RP\d+-|A[CL]\d{6}|CT[ADC]-|AP0\d)")
# biotypes worth a label when they have a real symbol; everything else is dropped unless
# it is protein-coding, which is always kept
NONCODING_OK = {"lncRNA", "antisense", "miRNA", "processed_transcript"}


def _get(url):
    with urllib.request.urlopen(url, timeout=60) as r:
        return json.load(r)


def fetch(host, chrom, lo, hi):
    return _get(f"{host}/overlap/region/human/{chrom.replace('chr', '')}:{lo}-{hi}"
                "?feature=gene;content-type=application/json")


def assert_build(host, assembly, rsid, bp, locus):
    """Refuse to write a track from the wrong assembly.

    A gene track is silently wrong when the build is wrong -- the genes are real and the
    coordinates are plausible, they just belong to another genome. So this is checked
    rather than assumed: the lead variant's position on `assembly` must equal the bp our
    own panel table carries for it."""
    m = [x for x in _get(f"{host}/variation/human/{rsid}?content-type=application/json")
         .get("mappings", []) if str(x.get("assembly_name", "")).lower() == assembly]
    if not m:
        raise SystemExit(f"{locus}: Ensembl returned no {assembly} mapping for {rsid}; "
                         "cannot confirm the build, refusing to write")
    got = int(m[0]["start"])
    if got != int(bp):
        raise SystemExit(
            f"{locus}: BUILD MISMATCH. Our table puts {rsid} at {int(bp):,}; Ensembl "
            f"{assembly} puts it at {got:,}. This pipeline is GRCh37 -- pass "
            "--assembly grch37 (or fix the panel table), do not write this track.")
    print(f"  build OK: {rsid} at {got:,} on {assembly}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--locus", nargs="+", required=True, help="chr1_238935959_239252559 ...")
    ap.add_argument("--locus-data", default=os.path.join(here, "..", "results", "fig2_loci"))
    ap.add_argument("--assembly", default="grch37", choices=sorted(HOST),
                    help="GRCh37 is correct for this pipeline; the check below enforces it")
    ap.add_argument("--flank-kb", type=float, default=300.0,
                    help="also keep named genes this far outside the interval, marked "
                         "flank=1, so a gene desert still gets its nearest neighbour")
    args = ap.parse_args()

    host = HOST[args.assembly]
    for L in args.locus:
        chrom, lo, hi = L.split("_")[0], int(L.split("_")[1]), int(L.split("_")[2])
        print(f"{L}:")
        r2 = pd.read_csv(os.path.join(args.locus_data, f"fig2_{L}_r2_to_lead.csv.gz"))
        lead = r2.lead_rsid.iloc[0]
        assert_build(host, args.assembly, lead,
                     r2.loc[r2.rsid == lead, "bp"].iloc[0], L)
        pad = int(args.flank_kb * 1000)
        g = [x for x in fetch(host, chrom, lo - pad, hi + pad)
             if x.get("external_name")
             and (x["biotype"] == "protein_coding"
                  or (x["biotype"] in NONCODING_OK
                      and not CLONE.match(x["external_name"])))]
        rows = [dict(gene=x["external_name"], start=x["start"], end=x["end"],
                     strand=x["strand"], biotype=x["biotype"],
                     flank=int(x["end"] < lo or x["start"] > hi)) for x in g]
        d = pd.DataFrame(rows).sort_values(["flank", "start"])
        inside = d[d.flank == 0]
        if (inside.biotype == "protein_coding").any():
            # the interval has protein-coding genes of its own; flanking ones are clutter
            d = inside
        else:
            # a coding desert -- chr1:238,935,959-239,252,559 holds two pseudogenes and
            # nothing else. Keep whatever is inside and add the nearest protein-coding
            # gene on each side, so the panel still says what the region is near.
            def pick(sub, key):
                c = sub[sub.biotype == "protein_coding"]
                if not len(c):
                    return c
                return c.loc[[c[key].idxmax() if key == "end" else c[key].idxmin()]]
            out_ = d[d.flank == 1]
            d = pd.concat([inside, pick(out_[out_.end < lo], "end"),
                           pick(out_[out_.start > hi], "start")])
        out = os.path.join(args.locus_data, f"fig2_{L}_genes.csv")
        d.to_csv(out, index=False)
        n_in = int((d.flank == 0).sum())
        print(f"  {len(d)} named genes ({n_in} inside, {len(d) - n_in} flanking) -> {out}")


if __name__ == "__main__":
    main()

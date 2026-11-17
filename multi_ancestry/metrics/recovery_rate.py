#!/usr/bin/env python3
"""Benchmark recovery rate by locus: does an arm's top variant match the published one?

This is the manuscript's Fig. 1b metric, rebuilt on the 2026-09 primary run and extended
from 5 arms to 9. It answers a different question from `plot_cs_split.py`: that one scores
every variant, this one scores every LOCUS, once.

The rule, unchanged from `manuscript_check/hierarchy_v3_gws_confident.py`:

  1. Stratum = the 122 loci where the published SuSiEx max PIP is >= 0.1 (`loci_confident122`),
     a subset of the 175 single-credible-set loci.
  2. An arm that does not resolve EXACTLY ONE credible set at a locus is counted as
     non-concordant. That is deliberately conservative and it is where most misses come from.
  3. Otherwise take the top-PIP member OF THAT ONE CREDIBLE SET -- not the locus maximum --
     and compare it to the published top variant.
  4. `exact` = same variant. `ld_match` = exact OR r2 > 0.95 to it. `ld_match` is the
     manuscript's headline; `exact` is reported alongside as the ordering check.

NO GWS FILTER, AND THAT IS NOW THE CORRECT CHOICE. The published
82.8/76.2/71.3/63.9/56.6% figures applied a "drop any credible set with no genome-wide
significant member" filter to OUR credible sets, as our reconstruction of Yuan et al.'s own
procedure. With their full output now in hand (`analysis/benchmark_susiex_yuan/`) that filter
does not apply to their credible sets, so applying it to ours would make the comparison
asymmetric -- penalising our arms with a gate the benchmark never passed through. Both sides
are therefore scored unfiltered. (Its three sumstats inputs were also lost with the old
workspace and Task B24 could not run it; that is no longer the reason it is absent.)

ONE DEPARTURE FROM THE PUBLISHED PIPELINE, deliberate and material:


**r2 comes from the full EUR panel, via `task_a_out/proxy_map_r095.csv.gz`.** The published
scorer looked r2 up in each arm's own harmonized npz and scored a miss when the published top
variant was absent from it. For FOUR loci the published lead marker is absent from our
reference LD entirely -- chr2_104919677, chr14_72357034, chr15_82425738, chr22_51078291 -- so
that lookup could not succeed there by construction, and those four were being charged to the
methods rather than to the panel. The proxy map resolves them through the full EUR panel. This
is the +4 that reconciled our recount with the manuscript, and it is why the Methods must say
the published lead marker is missing from our LD for these four loci. `--no-proxy-rescue`
reports the rate without it, so the size of the effect stays visible rather than assumed.

EMULATION: `--emulate-published` reproduces the manuscript's numbers so they can be traced.

It re-imposes the one thing the old scorer did that this one does not: `r2_rsid` and `r2_pos`
return None -- scored as a miss -- when EITHER variant is absent from the matrix r2 is looked
up in, so a proxy only counted when both our top variant and the published lead were in that
matrix. The flag applies exactly that gate, per arm, using the right matrix for each:

  * the three per-ancestry arms looked r2 up in the UNION harmonised npz (`NPZ_SOURCE` in
    hierarchy_v3_gws_confident.py points multisusie, susiex AND mesusie at the same file).
    That is MultiSuSiE's variant set, NOT SuSiEx's wider own set -- SuSiEx sees 11,051 more
    variants than the npz, so it was being gated by a matrix narrower than its own input.
    Missing this is why a first attempt left SuSiEx 4 loci high.
  * the merged arms looked r2 up in their own merged `.z` (`MATCH_BY_POS`).

Result, against the published counts out of 122:

  MultiSuSiE      101  ->  101   82.8%  EXACT
  MESuSiE          87  ->   87   71.3%  EXACT
  FINEMAP-merged   69  ->   69   56.6%  EXACT
  SuSiEx           94  vs   93   77.0%  one locus high
  SuSiE-merged     86  vs   78   70.5%  eight loci high

The last two are NOT reproducible from the committed data, for a reason worth stating plainly
rather than tuning away. The gate above is about the lookup matrix's VARIANT LIST; what
remains is the r2 VALUE it returns. The old scorer computed r2 inside the harmonised npz and
the merged `.z`; this scorer reads it from `proxy_map_r095.csv.gz`, built on the full EUR
panel. Merged LD is averaged across ancestries and attenuates pairwise r2, so eight of
SuSiE-merged's proxies clear 0.95 in the EUR panel and did not in the merged matrix. Closing
that would need the merged LD matrices themselves, which are not in the repo. SuSiEx's single
locus is the same effect at the harmonised npz, or the GWS filter; it cannot be separated
without the lost per-locus output.

So: three arms reproduce exactly, and the other two are explained by a mechanism that is
measured (8 and 1 loci) rather than assumed. The per-locus CSV carries `proxy` and
`proxy_blocked_by_lookup_domain` so any individual locus can be traced.


MERGED-LD r2 DOES NOT CLOSE THE GAP -- tested 2026-09-25, hypothesis NOT supported.

Task B32 supplied `r2_merged` for all 2,177 proxy pairs, and `--merged-r2` judges the
merged-LD arms by it instead of the EUR panel, which is what the published `r2_pos()` did.
The expectation was that this would reproduce SuSiE-merged's 63.9% exactly. It does not:

                        emulate (EUR r2)   emulate (merged r2)   published
  MultiSuSiE                   101                101              101
  MESuSiE                       87                 87               87
  SuSiEx                        94                 94               93
  SuSiE-merged                  86                 83               78
  FINEMAP-merged                69                 66               69

It closes 3 of SuSiE-merged's 8 and BREAKS FINEMAP-merged, which reproduced exactly before.
Only six pairs flip in total (three per arm), none of them already blocked by the
lookup-domain gate, so the effect is real but far too small and pointed the wrong way for
FINEMAP.

The leading explanation, NOT yet verified: Task B16 rebuilt the merged LD and the merged
FINEMAP inputs. The matrix exported here is therefore not necessarily the matrix the
manuscript's scorer read, and a partially-overlapping matrix is exactly what would close some
pairs, open others, and help one arm while hurting another. Confirming it needs the
pre-B16 merged LD, which may not have survived. Recorded as an open question rather than
tuned around: `--merged-r2` stays available and off by default, and the honest summary remains
three arms exact, two accounted for by mechanism but not reproduced.


RECONCILIATION WITH THE PUBLISHED FIGURES. The manuscript's rule was `ld_match = exact OR
r2 > 0.95`, applied at EVERY locus -- the proxy allowance was the matching rule itself, not a
fallback used only where the lead variant was missing. The defect is in WHERE r2 was looked
up: `hierarchy_v3_gws_confident.py` returns None when it cannot compute r2, and None is
scored as a miss, so "no proxy exists" and "could not check" are the same outcome.

Decomposing our counts against the published ones, on the assumption that exact matches agree
(they need no LD, and the runs are identical):

  arm             published = exact + proxy | ours = exact + proxy | extra proxies explained
  MultiSuSiE          101      95      6    |  105     95     10   | 4 = the four loci
  MESuSiE              87      83      4    |   91     83      8   | 4 = the four loci
  FINEMAP-merged       69      56     13    |   72     56     16   | 3 = the four loci
  SuSiE-merged         78      70      8    |   90     70     20   | 4 four loci + 8 merged-LD
  SuSiEx               93      89      4    |   97     89      8   | 2 four loci + 2 UNEXPLAINED

So three arms reproduce EXACTLY once the four loci are removed, and SuSiE-merged reproduces
exactly once the merged-LD attenuation below is also removed. SuSiEx leaves a residual of 2
loci. The remaining uncontrolled difference is the GWS filter, which the published pipeline
applied and we do not: it can only remove credible sets, so it can turn a match into a miss
(1 CS -> 0) or a forced miss into a scorable locus (2 CS -> 1), in either direction. It cannot
be rerun (its three sumstats inputs were lost with the old workspace, Task B24) and the
per-locus output that would let us diff locus by locus, `hierarchy_v3_gws_filtered_all175.csv`,
lived in the same lost CHECK_DIR. So those 2 loci are accounted for by mechanism but cannot be
named.


NO ARM'S RUN MOVED. Scoring the OLD (B12/B13/B16) exports and the B30 primary exports
through this scorer gives identical results for every arm tested -- SuSiEx 89 exact / 97
matched, MultiSuSiE 95/105, MESuSiE 83/91, SuSiE-merged 70/90 -- so the SuSiEx allele fix,
coverage 0.95 -> 0.99 and the rebuilt LD change none of these numbers. Every difference from
the published figures is the SCORER, and it decomposes exactly:

  MultiSuSiE    82.8% -> 86.1%   101 -> 105 loci   +4
  SuSiEx        76.2% -> 79.5%    93 ->  97 loci   +4
  MESuSiE       71.3% -> 74.6%    87 ->  91 loci   +4
  SuSiE-merged  63.9% -> 73.8%    78 ->  90 loci  +12
  FINEMAP-merged 56.6% -> 59.0%   69 ->  72 loci   +3

The three per-ancestry arms each gain EXACTLY the four loci whose published lead marker is
absent from our reference LD -- the +4 finding, and nothing else. SuSiE-merged gains those
four plus eight more, for the separate reason below.


THIS IS ALSO WHY SuSiE-merged MOVED 63.9% -> 73.8%, AND IT IS NOT THE RUN. Scoring B13's old
run, B16's rebuilt-LD run and B30's primary run through this scorer gives identical results --
13/108/1 credible-set counts, 70 exact, 90 matched, 73.8% -- so neither the SuSiEx allele fix
(which touches only SuSiEx's own join), nor coverage 0.95 -> 0.99, nor the rebuilt LD moves
this arm at all. The entire gain is the r2 lookup: published 70 exact + 8 proxy rescues = 78,
here 70 + 20 = 90.

The mechanism is the r2 VALUE, not a narrow matrix. The published scorer routed the merged
arms through `MATCH_BY_POS` -> `r2_pos()`, which computes r2 in the locus's own MERGED
`.z`/`.ld`, while the per-ancestry arms used `r2_rsid()` against the harmonised union npz's
EUR slice. It is tempting to explain the gap by the merged matrix being narrower -- but it is
not: the merged arms carry a median 405 variants per locus against MultiSuSiE's 408, and at
16 of the 20 rescued loci the published lead IS in SuSiE-merged's own domain, so the old
lookup could have run and still did not clear 0.95. What differs is the number it returns.
Merged LD is averaged across ancestries, which attenuates pairwise r2 relative to a
single-ancestry panel, so a proxy pair that clears 0.95 in the EUR panel need not clear it in
the merged matrix. The published scoring therefore held the merged arms to a stricter proxy
test than the per-ancestry arms, using the same nominal threshold. Only the remaining 4 of
those 20 are domain failures, and they are the four loci named above.

THE DOMAIN CEILING IS ITS OWN RESULT, and it is MESuSiE's, not merged LD's. `lead_in_domain`
counts the loci where the published lead variant is present in the arm's input at all; an arm
cannot exact-match a variant it never saw. SuSiEx and the EUR-only arms see it at 121/122
(99.2%), MultiSuSiE and every merged arm at 116/122 (95.1%), but MESuSiE -- which fine-maps
the EUR/EAS INTERSECTION rather than the union -- at only 109/122 (89.3%), off a median 281
variants per locus against everyone else's 405-426. So roughly a third of MESuSiE's apparent
shortfall against MultiSuSiE is a ceiling imposed by its input, before its model is consulted.
Merged LD does NOT share this shortcoming; the intersection is what causes it.

Usage:
  python3 recovery_rate.py --export ../aou_ld/task_b30_out --out figures/recovery_122.csv
"""
import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bench_sensfdr import BENCH, bench_cs_loci, KAI

HERE = os.path.dirname(os.path.abspath(__file__))
PROXY = os.path.join(HERE, "..", "aou_ld", "task_a_out", "proxy_map_r095.csv.gz")
PROXY_MERGED = os.path.join(HERE, "..", "aou_ld", "task_b32_out",
                            "proxy_map_r095_with_merged.csv.gz")
# the published scorer sent these two down r2_pos(), computing r2 inside the merged matrix
MERGED_ARMS = {"SuSiE-merged", "FINEMAP-merged", "RSparsePro", "CARMA"}
# "193" is not a file: it is every locus where the benchmark formed at least one credible
# set, the same stratum the variant-level figures use. Needed because the 13 variants with
# benchmark PIP > 0.9 do not all fall inside the 122 (one locus has >=2 benchmark CSs).
STRATA = {"122": "loci_confident122.txt", "175": "loci_single175.txt", "193": None}
# The grouping, same order as the CS-split figure so the two read together.
# SuSiEx-noambig (Task B33) is the ambiguous-free rerun; it sits next to SuSiEx-ours so
# the two can be read against each other, and is simply absent when B33 has not been
# exported.
ORDER = ["SuSiEx-ours", "SuSiEx-noambig", "MultiSuSiE", "MESuSiE", "FINEMAP-merged",
         "EUR-FINEMAP",
         "SuSiE-merged", "EUR-SuSiE-RSS", "RSparsePro", "CARMA"]
# Published lead marker absent from our reference LD; see the docstring.
LEAD_ABSENT = ["chr2_104919677", "chr14_72357034", "chr15_82425738", "chr22_51078291"]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--export", default=os.path.join(HERE, "..", "aou_ld", "task_v2_out", "noambig"))
    ap.add_argument("--stratum", default="122", choices=sorted(STRATA))
    ap.add_argument("--purity", type=float, default=None,
                    help="gate credible sets on common_purity >= this before counting them; "
                         "omit to use every credible set the arm reported")
    ap.add_argument("--min-pip", type=float, default=None,
                    help="also require the ARM's own PIP at its top variant to exceed this, "
                         "so both sides of the comparison must be confidently resolved "
                         "(the 122-locus stratum already imposes it on the benchmark side)")
    ap.add_argument("--emulate-published", action="store_true",
                    help="reproduce the manuscript's numbers by re-imposing the old scorer's "
                         "lookup-domain gate: a proxy counts only if BOTH our top variant "
                         "and the published lead are in the matrix r2 was looked up in. "
                         "See the EMULATION section of the docstring for what it does and "
                         "does not reproduce.")
    ap.add_argument("--drop-palindromic", action="store_true",
                    help="drop A/T and C/G variants from OUR side before scoring, so we are "
                         "judged on the same variant domain the benchmark had. Yuan et al.'s SuSiEx "
                         "input was a 'SumStats.noAmbig' directory -- ambiguous SNPs were "
                         "stripped upstream -- and our set minus palindromic reproduces his "
                         "EXACTLY at all 193 loci (146,698 = 146,698). Without this an arm "
                         "whose top CS variant is palindromic is a forced miss.")
    ap.add_argument("--merged-r2", action="store_true",
                    help="for the merged-LD arms, judge proxies by r2 computed in the MERGED "
                         "matrix (task_b32_out) instead of the EUR panel, which is what the "
                         "published scorer's r2_pos() did. Combine with --emulate-published.")
    ap.add_argument("--proxy-only-if-lead-absent", action="store_true",
                    help="allow a proxy ONLY where the published lead is absent from the "
                         "arm's input. Where the lead IS present, the arm saw it, weighed it "
                         "and preferred something else -- scoring that as recovery credits "
                         "the method for a choice it did not make. This is the correction "
                         "that stops merged LD being inflated.")
    ap.add_argument("--no-proxy-rescue", action="store_true",
                    help="score exact matches only, no r2>0.95 proxies")
    ap.add_argument("--suffix", default="-noambig", help="export label suffix, e.g. --suffix=-noambig for task_v2_out/noambig")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    if STRATA[args.stratum] is None:
        loci = sorted(bench_cs_loci(KAI)[0])
    else:
        loci = sorted({l.strip() for l in
                       open(os.path.join(HERE, STRATA[args.stratum])) if l.strip()})

    # --- the published top variant per locus -----------------------------------------
    b = pd.read_csv(BENCH)
    b = b[b.locus.isin(loci)]
    top = b.loc[b.groupby("locus").pip.idxmax()][["locus", "rsid", "bp", "pip"]]
    top = top.rename(columns={"rsid": "pub_rsid", "bp": "pub_bp", "pip": "pub_pip"})
    top = top.set_index("locus")
    missing = [l for l in loci if l not in top.index]
    if missing:
        raise SystemExit(f"[error] no benchmark PIPs for {len(missing)} loci: {missing[:3]}")

    # --- proxies of that variant, from the full EUR panel -----------------------------
    pm = pd.read_csv(PROXY_MERGED if args.merged_r2 else PROXY)
    prox, prox_merged = {}, {}
    for l, g in pm[pm.locus.isin(loci)].groupby("locus"):
        lead = top.at[l, "pub_rsid"]
        gg = g[(g.rsid_a == lead) & (g.r2 > 0.95)]
        # the map is emitted once per unordered pair, so accept the lead on either side
        gg2 = g[(g.rsid_b == lead) & (g.r2 > 0.95)]
        prox[l] = (set(gg.bp_b.astype(int)) | set(gg2.bp_a.astype(int)))
        if args.merged_r2:
            m1 = g[(g.rsid_a == lead) & (g.r2_merged > 0.95)]
            m2 = g[(g.rsid_b == lead) & (g.r2_merged > 0.95)]
            prox_merged[l] = set(m1.bp_b.astype(int)) | set(m2.bp_a.astype(int))

    rows, per_locus = [], []
    for arm in ORDER:
        f = os.path.join(args.export, f"primary_{arm}{args.suffix}_cs_purity.csv.gz")
        if not os.path.exists(f):
            print(f"[warn] no export for {arm}, skipped")
            continue
        d = pd.read_csv(f, low_memory=False)
        # The arm's own variant DOMAIN, from its full pip table: a locus whose published
        # lead variant is not in the arm's input can never be an exact match, whatever the
        # method does. This is a ceiling imposed by the input, not by the fine-mapper, and
        # it differs sharply by arm -- MESuSiE takes the EUR/EAS INTERSECTION and is ~30%
        # narrower than everything else.
        # The old scorer's r2 lookup matrix, per arm: NPZ_SOURCE points the three
        # per-ancestry arms at the UNION harmonised npz -- which is MultiSuSiE's variant
        # set, NOT SuSiEx's wider own set -- and MATCH_BY_POS points the merged arms at
        # their own merged .z. That distinction is what makes SuSiEx's emulation work.
        lookup = None
        if args.emulate_published:
            src = ("MultiSuSiE" if arm in ("SuSiEx-ours", "MultiSuSiE", "MESuSiE") else arm)
            lf = os.path.join(args.export, f"primary_{src}_pip.csv.gz")
            if os.path.exists(lf):
                lv = pd.read_csv(lf, low_memory=False)
                lv = lv[lv.locus.isin(loci)]
                lookup = {l: set(g.bp.astype(int)) for l, g in lv.groupby("locus")}
        pf = os.path.join(args.export, f"primary_{arm}{args.suffix}_pip.csv.gz")
        dom_hit = np.nan
        if os.path.exists(pf):
            pv = pd.read_csv(pf, low_memory=False)
            pv = pv[pv.locus.isin(loci)]
            dom = {l: set(g.bp.astype(int)) for l, g in pv.groupby("locus")}
            dom_hit = sum(int(top.at[l, "pub_bp"]) in dom.get(l, set()) for l in loci)
            med_n = int(pv.groupby("locus").size().median())
        d = d[d.locus.isin(loci) & d.cs_ids_all.notna()].copy()
        if args.drop_palindromic:
            PAL = [frozenset(("A", "T")), frozenset(("C", "G"))]
            keep_row = [not (pd.notna(x) and pd.notna(y)
                             and frozenset((str(x).upper(), str(y).upper())) in PAL)
                        for x, y in zip(d.a1, d.a2)]
            d = d[keep_row]
        if args.purity is not None:
            d = d[d.common_purity.notna() & (d.common_purity >= args.purity)]
        # a variant can belong to more than one credible set
        d["cs"] = d.cs_ids_all.astype(str).str.split(r"[;,|]")
        d = d.explode("cs")
        d["cs"] = d.cs.str.strip()

        n_matched = n_exact = n_exact_hi = n_exact_lo = 0
        n_cs_counts = {0: 0, 1: 0, 2: 0}
        for l in loci:
            g = d[d.locus == l]
            ncs = g.cs.nunique()
            n_cs_counts[0 if ncs == 0 else (1 if ncs == 1 else 2)] += 1
            rec = dict(arm=arm, locus=l, n_cs=ncs, pub_rsid=top.at[l, "pub_rsid"],
                       pub_bp=int(top.at[l, "pub_bp"]), top_rsid=None, top_bp=None,
                       exact=False, ld_match=False, proxy=False,
                       proxy_blocked_by_lookup_domain=False)
            if ncs == 1:
                r = g.loc[g.pip.idxmax()]
                # A confident call on the arm's side too, when --min-pip is given: the
                # stratum already requires the BENCHMARK's max PIP >= 0.1, so without this
                # the test is one-sided -- an arm could "match" on a top variant it holds
                # at PIP 0.02 and be scored the same as one holding it at 0.99.
                confident = args.min_pip is None or float(r.pip) > args.min_pip
                is_lead = int(r.bp) == int(top.at[l, "pub_bp"])
                # split exact matches by OUR confidence in the lead, independent of
                # --min-pip, so the cost of demanding confidence is visible either way
                if is_lead:
                    if float(r.pip) >= 0.1:
                        n_exact_hi += 1
                    else:
                        n_exact_lo += 1
                exact = confident and is_lead
                pset = (prox_merged if (args.merged_r2 and arm in MERGED_ARMS)
                        else prox).get(l, set())
                is_prox = int(r.bp) in pset
                if args.proxy_only_if_lead_absent and dom is not None:
                    # the lead is in our input: we saw it and preferred another variant
                    if int(top.at[l, "pub_bp"]) in dom.get(l, set()):
                        is_prox = False
                blocked = False
                if is_prox and lookup is not None:
                    dm = lookup.get(l, set())
                    # r2_rsid/r2_pos return None -- scored as a miss -- when EITHER variant
                    # is absent from the lookup matrix, so both must be present
                    blocked = not (int(top.at[l, "pub_bp"]) in dm and int(r.bp) in dm)
                ld = exact or (confident and not args.no_proxy_rescue and is_prox
                               and not blocked)
                rec["proxy"] = bool(is_prox and not exact)
                rec["proxy_blocked_by_lookup_domain"] = bool(blocked)
                rec.update(top_rsid=r.rsid, top_bp=int(r.bp), exact=bool(exact),
                           ld_match=bool(ld))
                n_exact += bool(exact)
                n_matched += bool(ld)
            per_locus.append(rec)
        n = len(loci)
        rows.append(dict(arm=arm, n_loci=n, n_cs_0=n_cs_counts[0], n_cs_1=n_cs_counts[1],
                         n_cs_ge2=n_cs_counts[2], exact=n_exact, ld_match=n_matched,
                         exact_rate=n_exact / n, recovery_rate=n_matched / n,
                         median_variants=med_n, lead_in_domain=dom_hit,
                         domain_ceiling=dom_hit / n,
                         # exact lead matches split by OUR PIP at that variant
                         exact_lead_hi=n_exact_hi, exact_lead_lo=n_exact_lo,
                         # conditional denominator: loci where the lead is in our input at
                         # all. An arm cannot exact-match a variant it never saw, so this
                         # separates the model's failure from its input's.
                         rate_lead_present=((n_exact_hi + n_exact_lo) / dom_hit
                                            if dom_hit else np.nan)))
    R = pd.DataFrame(rows)
    P = pd.DataFrame(per_locus)

    gate = "none" if args.purity is None else f"common_purity >= {args.purity}"
    mp = "none" if args.min_pip is None else f"arm PIP > {args.min_pip:g}"
    print(f"\nStratum {args.stratum} ({len(loci)} loci) | purity gate: {gate} | "
          f"arm confidence: {mp} | "
          f"proxy rescue: {'OFF' if args.no_proxy_rescue else 'ON (r2>0.95, full EUR panel)'}")
    print("GWS filter not applied, to either side (see the docstring)\n")
    show = R.copy()
    show["recovery"] = show.recovery_rate.map(lambda v: f"{v:.1%}")
    show["exact"] = show.exact_rate.map(lambda v: f"{v:.1%}")
    show["ceiling"] = show.domain_ceiling.map(lambda v: f"{v:.1%}")
    show["cond"] = show.rate_lead_present.map(lambda v: f"{v:.1%}")
    print(show[["arm", "median_variants", "lead_in_domain", "ceiling", "n_cs_0", "n_cs_1",
                "n_cs_ge2", "exact_lead_hi", "exact_lead_lo", "ld_match", "recovery",
                "exact", "cond"]].to_string(index=False))

    sub = P[P.locus.str.startswith(tuple(x + "_" for x in LEAD_ABSENT))]
    if len(sub):
        print(f"\nThe {sub.locus.nunique()} loci whose published lead marker is absent from "
              f"our reference LD:")
        print(sub[sub.n_cs == 1].groupby("locus").agg(
            arms_scored=("arm", "size"), matched=("ld_match", "sum"),
            exact=("exact", "sum")).to_string())

    # --- Cochran's Q over the arms, on the SAME 122 loci (a within-locus repeated-measures
    # design, which is the point: the arms are scored on identical loci, so the paired
    # structure must be respected rather than treating nine rates as independent).
    from scipy.stats import chi2
    for label, sel in (("all 9 arms", ORDER),
                       ("the manuscript's 5", ["MultiSuSiE", "SuSiEx-ours", "MESuSiE",
                                               "SuSiE-merged", "FINEMAP-merged"])):
        sel = [a for a in sel if a in set(P.arm)]
        M = (P[P.arm.isin(sel)].pivot(index="locus", columns="arm", values="ld_match")
             .astype(int)[sel])
        k = M.shape[1]
        Lsum, Gsum = M.sum(axis=1).to_numpy(), M.sum(axis=0).to_numpy()
        denom = k * Lsum.sum() - (Lsum ** 2).sum()
        if denom > 0:
            Q = (k - 1) * (k * (Gsum ** 2).sum() - Gsum.sum() ** 2) / denom
            print(f"Cochran's Q ({label}, k={k}, n=122): Q={Q:.1f}, df={k-1}, "
                  f"P={chi2.sf(Q, k - 1):.2g}")

    if args.out:
        R.to_csv(args.out, index=False)
        P.to_csv(args.out.rsplit(".", 1)[0] + "_per_locus.csv", index=False)
        print(f"\n[done] wrote {args.out} and its _per_locus.csv")


if __name__ == "__main__":
    main()

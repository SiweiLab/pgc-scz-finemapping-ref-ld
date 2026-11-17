#!/bin/bash
# B52 §1 driver: harmonize_v2.py once per locus under /usr/bin/time -v (wall + peak RSS per locus),
# 6 in parallel, chr3_48182866 first. A failing locus exits 255 so xargs stops launching new ones.
# Usage: run_harmonize_b52.sh primary|noambig [extra harmonize_v2 args, e.g. --arms multisusie]  (B54 added the extra args)
set -u
V=$1; shift; EXTRA="$*"; OUT=~/work/v2_harmonized/$V; LOG=$OUT/_time_${TAG:-b52}; mkdir -p $LOG
FLAG=""; [ "$V" = noambig ] && FLAG="--drop-ambiguous"
export V2_LD_ROOT=~/work/task_b50/ld_matrices V2_SS_ROOT=~/work/task_b50/PGC_SCZ_sumstats OUT LOG FLAG EXTRA
cd "$(dirname "$0")"
( grep chr3_48182866 v2_inputs/loci_250.txt; grep -v chr3_48182866 v2_inputs/loci_250.txt ) |
  xargs -P 6 -I{} bash -c '/usr/bin/time -v -o $LOG/{}.time python harmonize_v2.py --loci {} --out $OUT $FLAG $EXTRA > $LOG/{}.out 2>&1 || { echo "FAIL {}" >> $LOG/FAILED; exit 255; }'
echo "exit $?"

#!/bin/bash
# Driver: harmonize.py once per locus under /usr/bin/time -v (wall + peak RSS per locus),
# 6 in parallel, chr3_48182866 first. A failing locus exits 255 so xargs stops launching new ones.
# Usage: run_harmonize.sh primary|noambig [extra harmonize args, e.g. --arms multisusie]
set -u
V=$1; shift; EXTRA="$*"; OUT=~/work/v2_harmonized/$V; LOG=$OUT/_time_${TAG:-run}; mkdir -p $LOG
FLAG=""; [ "$V" = noambig ] && FLAG="--drop-ambiguous"
export LD_ROOT=~/work/ld_matrices SS_ROOT=~/work/PGC_SCZ_sumstats OUT LOG FLAG EXTRA
cd "$(dirname "$0")"
( grep chr3_48182866 ../inputs/loci_250.txt; grep -v chr3_48182866 ../inputs/loci_250.txt ) |
  xargs -P 6 -I{} bash -c '/usr/bin/time -v -o $LOG/{}.time python harmonize.py --loci {} --out $OUT $FLAG $EXTRA > $LOG/{}.out 2>&1 || { echo "FAIL {}" >> $LOG/FAILED; exit 255; }'
echo "exit $?"

#!/bin/bash
# Unattended overnight finish of the CARMA L sweep. Launch DETACHED once the six CARMA sweep drivers are running:
#   setsid nohup bash overnight_lsweep_carma.sh &> ~/work/overnight_lsweep_carma.log &
#   1. wait until every CARMA sweep driver (run_finemap.py --arms CARMA --L ...) has exited
#   2. export each (variant, L) -> results/finemap/L_sweep[_with_ambiguous]/primary_CARMA[-noambig]_L{n}_*.csv.gz
#   3. write a machine summary
#   4. check the bucket backups (l_sweep_v2_2026-10/<V>/CARMA_L{n}) against local
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"; cd "$HERE"
R="$HERE/../results/finemap"
LOG=~/work/overnight_lsweep_carma; mkdir -p "$LOG"
say() { echo "[$(date -u +%FT%TZ)] $*" | tee -a "$LOG/steps.log"; }
B=${WORKSPACE_BUCKET:?set WORKSPACE_BUCKET to your workspace bucket}/manuscript_final_2026-09/fine_mapping_runs/l_sweep_v2_2026-10
: "${LD_ROOT:?set LD_ROOT to the local ld_matrices (ld_per_anc_eur_gctb + ld_per_anc_eas_aou) before launching}"

say "1. waiting for the CARMA sweep drivers"
while pgrep -f "run_finemap.py --arms CARMA --L" > /dev/null; do sleep 300; done
for V in noambig primary; do for L in 1 3 5; do
  n=$(ls ~/work/primary_v2_$V/CARMA_L$L/*_carma_pip.tsv 2>/dev/null | wc -l)
  say "$V L=$L: $n/250 loci with output"
done; done

say "2. export"
for V in noambig primary; do for L in 1 3 5; do
  HARMONIZED_INPUTS=~/work/v2_harmonized/$V python3 export_results.py --L $L --arms CARMA &> "$LOG/export_${V}_L$L.log" \
    && say "export $V L=$L ok" || say "export $V L=$L FAILED"
done; done
python3 - "$R" > "$R/L_sweep/overnight_lsweep_carma_summary.txt" <<'EOF'
import glob, sys
import pandas as pd
R = sys.argv[1]
print("file\trows\tloci\tpip>0.5\tpip>0.95")
for f in sorted(glob.glob(f"{R}/L_sweep*/primary_CARMA*_pip.csv.gz")):
    d = pd.read_csv(f, usecols=["locus", "pip"])
    print(f"{f}\t{len(d)}\t{d.locus.nunique()}\t{int((d.pip > 0.5).sum())}\t{int((d.pip > 0.95).sum())}")
EOF
cp "$LOG"/steps.log "$R/L_sweep/overnight_lsweep_carma_steps.log"

# (As run, the script committed these exports to the working repository here; omitted in this release.)

say "3. backups"
for V in noambig primary; do for L in 1 3 5; do
  src=~/work/primary_v2_$V/CARMA_L$L/
  gsutil -m -q rsync -r -c -x '(^|.*/)_inputs/.*' "$src" "$B/$V/CARMA_L$L/" >> "$LOG/backup.log" 2>&1
  ops=$(gsutil -m rsync -r -n -c -x '(^|.*/)_inputs/.*' "$src" "$B/$V/CARMA_L$L/" 2>&1 | grep -c "Would copy")
  say "backup $V L=$L: $ops differences after sync"
done; done
cp "$LOG"/steps.log "$R/L_sweep/overnight_lsweep_carma_steps.log"
# (As run, the script committed these exports to the working repository here; omitted in this release.)
# (As run, the script then stopped the cluster; omitted in this release.)

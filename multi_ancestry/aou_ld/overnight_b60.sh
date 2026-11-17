#!/bin/bash
# B60, unattended overnight finish. Launch DETACHED once the six CARMA sweep drivers are running:
#   setsid nohup bash overnight_b60.sh &> ~/work/overnight_b60.log &
#   1. wait until every CARMA sweep driver (run_primary_2026.py --arms CARMA --L ...) has exited
#   2. export each (variant, L) -> task_v2_out/lsweep/<V>/primary_CARMA[-noambig]_L{n}_*.csv.gz
#   3. write a machine summary, commit, push
#   4. check the bucket backups (l_sweep_v2_2026-10/<V>/CARMA_L{n}) against local, then STOP (not delete) this cluster
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"; cd "$HERE"
REPO="$(git rev-parse --show-toplevel)"
LOG=~/work/overnight_b60; mkdir -p "$LOG"
say() { echo "[$(date -u +%FT%TZ)] $*" | tee -a "$LOG/steps.log"; }
B=${WORKSPACE_BUCKET:?set WORKSPACE_BUCKET to your workspace bucket}/manuscript_final_2026-09/fine_mapping_runs/l_sweep_v2_2026-10
: "${V2_LD_ROOT:?set V2_LD_ROOT to the local ld_matrices (ld_per_anc_eur_gctb + ld_per_anc_eas_aou) before launching}"

say "1. waiting for the CARMA sweep drivers"
while pgrep -f "run_primary_2026.py --arms CARMA --L" > /dev/null; do sleep 300; done
for V in noambig primary; do for L in 1 3 5; do
  n=$(ls ~/work/primary_v2_$V/CARMA_L$L/*_carma_pip.tsv 2>/dev/null | wc -l)
  say "$V L=$L: $n/250 loci with output"
done; done

say "2. export"
for V in noambig primary; do for L in 1 3 5; do
  PRIMARY_V2_INPUTS=~/work/v2_harmonized/$V python3 export_primary_v2.py --L $L --arms CARMA &> "$LOG/export_${V}_L$L.log" \
    && say "export $V L=$L ok" || say "export $V L=$L FAILED"
done; done
python3 - > task_v2_out/lsweep/overnight_b60_summary.txt <<'EOF'
import glob
import pandas as pd
print("file\trows\tloci\tpip>0.5\tpip>0.95")
for f in sorted(glob.glob("task_v2_out/lsweep/*/primary_CARMA*_pip.csv.gz")):
    d = pd.read_csv(f, usecols=["locus", "pip"])
    print(f"{f}\t{len(d)}\t{d.locus.nunique()}\t{int((d.pip > 0.5).sum())}\t{int((d.pip > 0.95).sum())}")
EOF
cp "$LOG"/steps.log task_v2_out/lsweep/overnight_b60_steps.log
for V in noambig primary; do for L in 1 3 5; do cp ~/work/primary_v2_$V/state_L$L.json task_v2_out/lsweep/$V/ 2>/dev/null; done; done

say "3. commit, push"
cd "$REPO"
git add -f Post-hoc_analyses/fine-mapping/multi_ancestry/aou_ld/task_v2_out/lsweep   # -f: exports are *.csv.gz, ignored by the global *.* rule
git commit -q -m "B60 overnight (unattended): CARMA L sweep 1/3/5 seeded, both variants; exports + summary in task_v2_out/lsweep/" || say "nothing to commit"
for i in 1 2 3; do git pull --no-rebase -q origin fine-mapping && git push -q origin fine-mapping && { say "pushed $(git rev-parse --short HEAD)"; break; }; sleep 60; done

say "4. backups"
for V in noambig primary; do for L in 1 3 5; do
  src=~/work/primary_v2_$V/CARMA_L$L/
  gsutil -m -q rsync -r -c -x '(^|.*/)_inputs/.*' "$src" "$B/$V/CARMA_L$L/" >> "$LOG/backup.log" 2>&1
  ops=$(gsutil -m rsync -r -n -c -x '(^|.*/)_inputs/.*' "$src" "$B/$V/CARMA_L$L/" 2>&1 | grep -c "Would copy")
  say "backup $V L=$L: $ops differences after sync"
done; done
cp "$LOG"/steps.log Post-hoc_analyses/fine-mapping/multi_ancestry/aou_ld/task_v2_out/lsweep/overnight_b60_steps.log
git add -f Post-hoc_analyses/fine-mapping/multi_ancestry/aou_ld/task_v2_out/lsweep/overnight_b60_steps.log
git commit -q -m "B60 overnight: backup check log" && git pull --no-rebase -q origin fine-mapping && git push -q origin fine-mapping
NAME=$(curl -s -H "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/instance/attributes/dataproc-cluster-name)
REGION=$(curl -s -H "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/instance/attributes/dataproc-region)
say "stopping cluster $NAME ($REGION)"
gcloud dataproc clusters stop "$NAME" --region "$REGION" --async >> "$LOG/steps.log" 2>&1 || say "stop FAILED: stop it by hand"

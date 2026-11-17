#!/bin/bash
# B61, unattended overnight chain. Launch DETACHED:  setsid nohup bash overnight_b61.sh &> ~/work/overnight_b61.log &
#   1. wait for both seeded chr3 CARMA fits (B59) to exit
#   2. export the final v2 primary run, all 9 arms, both variants  -> task_v2_out/{noambig,primary}/   (gate 7 inputs)
#   3. L sweep of the six non-CARMA sweep arms, both variants (CARMA's sweep runs separately, B60)
#      L = 1, 3: SuSiEx MultiSuSiE MESuSiE SuSiE-merged FINEMAP-merged RSparsePro
#      L = 5:    the same minus FINEMAP-merged (its primary IS k = 5)
#      L = 10:   SuSiEx (the sweep runs mult-step off, a different fit) and FINEMAP-merged
#   4. export the sweep -> task_v2_out/lsweep/<V>/, write a machine summary, commit, push
#   5. check the bucket backups exist, then STOP (not delete) this cluster
# Every setting comes from run_params through run_primary_2026.py; nothing is changed per locus. Failures are logged, never retried.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"; cd "$HERE"
REPO="$(git rev-parse --show-toplevel)"
LOG=~/work/overnight_b61; mkdir -p "$LOG"
say() { echo "[$(date -u +%FT%TZ)] $*" | tee -a "$LOG/steps.log"; }
export V2_LD_ROOT="${V2_LD_ROOT:-$HOME/work/task_b50/ld_matrices}"
B=${WORKSPACE_BUCKET:?set WORKSPACE_BUCKET to your workspace bucket}/manuscript_final_2026-09/fine_mapping_runs/l_sweep_v2_2026-10
JOBS=7   # per driver; two drivers (one per variant) at a time on the n2-standard-16 master

say "1. waiting for chr3 CARMA"
while pgrep -f "run_carma_pilot_b29.R.*chr3_48182866_53553354" > /dev/null || pgrep -f "run_primary_2026.py.*--loci chr3_48182866_53553354" > /dev/null; do sleep 300; done
for V in noambig primary; do
  f=~/work/primary_v2_$V/CARMA/chr3_48182866_53553354_carma_pip.tsv
  [ -s "$f" ] && say "chr3 $V CARMA done: $(($(wc -l < $f)-1)) rows" || say "WARNING chr3 $V CARMA has NO output"
done

say "2. final export, 9 arms x 2"
for V in noambig primary; do
  PRIMARY_V2_INPUTS=~/work/v2_harmonized/$V python3 export_primary_v2.py &> "$LOG/export_$V.log" && say "export $V ok" || say "export $V FAILED"
done

say "3. L sweep, non-CARMA arms"
run_L() {  # $1 = L, $2 = driver arms
  for V in noambig primary; do
    PRIMARY_V2_INPUTS=~/work/v2_harmonized/$V python3 run_primary_2026.py --arms "$2" --L "$1" --jobs $JOBS \
        --backup-prefix "$B/$V" &> "$LOG/sweep_${V}_L$1.log" &
  done
  wait
  say "sweep L=$1 done: $(grep -h -E '^\s+\S+\s+ok=' "$LOG"/sweep_*_L$1.log | tr -s ' ' | paste -sd ';')"
}
run_L 1  "SuSiEx,MultiSuSiE,MESuSiE,SuSiE-merged,FINEMAP-merged,RSparsePro"
run_L 3  "SuSiEx,MultiSuSiE,MESuSiE,SuSiE-merged,FINEMAP-merged,RSparsePro"
run_L 5  "SuSiEx,MultiSuSiE,MESuSiE,SuSiE-merged,RSparsePro"
run_L 10 "SuSiEx,FINEMAP-merged"

say "4. export the sweep"
exp_L() {  # $1 = L, $2 = export labels
  for V in noambig primary; do
    PRIMARY_V2_INPUTS=~/work/v2_harmonized/$V python3 export_primary_v2.py --L "$1" --arms "$2" &> "$LOG/export_${V}_L$1.log" \
      && say "export $V L=$1 ok" || say "export $V L=$1 FAILED"
  done
}
exp_L 1  "SuSiEx-ours,MultiSuSiE,MESuSiE,SuSiE-merged,FINEMAP-merged,RSparsePro"
exp_L 3  "SuSiEx-ours,MultiSuSiE,MESuSiE,SuSiE-merged,FINEMAP-merged,RSparsePro"
exp_L 5  "SuSiEx-ours,MultiSuSiE,MESuSiE,SuSiE-merged,RSparsePro"
exp_L 10 "SuSiEx-ours,FINEMAP-merged"

python3 - > task_v2_out/overnight_b61_summary.txt <<'EOF'
import glob, gzip, os
import pandas as pd
print("file\trows\tloci\tpip>0.5\tpip>0.95")
for f in sorted(glob.glob("task_v2_out/*/primary_*_pip.csv.gz") + glob.glob("task_v2_out/lsweep/*/primary_*_pip.csv.gz")):
    d = pd.read_csv(f, usecols=["locus", "pip"])
    print(f"{f}\t{len(d)}\t{d.locus.nunique()}\t{int((d.pip > 0.5).sum())}\t{int((d.pip > 0.95).sum())}")
EOF
cp "$LOG"/steps.log task_v2_out/overnight_b61_steps.log

say "5. commit, push, check backups, stop"
cd "$REPO"
# [as run; disabled in this release: committed results to the working repo / stopped the cluster] git add -f Post-hoc_analyses/fine-mapping/multi_ancestry/aou_ld/task_v2_out   # -f: exports are *.csv.gz, ignored by the global *.* rule
# [as run; disabled in this release: committed results to the working repo / stopped the cluster] git commit -q -m "B61 overnight (unattended): final v2 exports both variants (seeded CARMA incl. chr3) + L sweep of the six non-CARMA arms; see task_v2_out/overnight_b61_steps.log, _summary.txt" || say "nothing to commit"
# [as run; disabled in this release: committed results to the working repo / stopped the cluster] for i in 1 2 3; do git pull --no-rebase -q origin fine-mapping && git push -q origin fine-mapping && { say "pushed $(git rev-parse --short HEAD)"; break; }; sleep 60; done
for V in noambig primary; do say "bucket $V: $(gsutil ls "$B/$V/" 2>/dev/null | wc -l) arm dirs under $B/$V/"; done
say "primary backup CARMA $(gsutil ls ${WORKSPACE_BUCKET:?set WORKSPACE_BUCKET to your workspace bucket}/manuscript_final_2026-09/fine_mapping_runs/primary_v2_2026-09/noambig/CARMA/ 2>/dev/null | grep -c _carma_pip.tsv) noambig pip files"
# [as run; disabled in this release: committed results to the working repo / stopped the cluster] NAME=$(curl -s -H "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/instance/attributes/dataproc-cluster-name)
# [as run; disabled in this release: committed results to the working repo / stopped the cluster] REGION=$(curl -s -H "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/instance/attributes/dataproc-region)
# [as run; disabled in this release: committed results to the working repo / stopped the cluster] say "stopping cluster $NAME ($REGION)"
# [as run; disabled in this release: committed results to the working repo / stopped the cluster] gcloud dataproc clusters stop "$NAME" --region "$REGION" --async >> "$LOG/steps.log" 2>&1 || say "stop FAILED: stop it by hand"

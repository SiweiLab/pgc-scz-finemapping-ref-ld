#!/usr/bin/env bash
# Build the fine-mapping env from finemap_env.yml, then install the pieces conda
# cannot hold. Every version is pinned to what the runs used (VERSIONS.md). Idempotent.
#
# Usage: bash environment/setup_env.sh [prefix]      # default prefix ~/envs/finemap_env
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PREFIX="${1:-$HOME/envs/finemap_env}"
TOOLS="${TOOLS_DIR:-$HOME/tools}"

# 1. conda env at a USER-WRITABLE prefix. The Dataproc base env is root-owned, which is
#    what blocked an earlier CARMA install (EnvironmentNotWritableError).
if [ ! -x "$PREFIX/bin/Rscript" ]; then
  mamba env create -y -p "$PREFIX" -f "$HERE/finemap_env.yml"
fi
R="$PREFIX/bin/Rscript"

# 2. GitHub-only R packages at pinned commits
"$R" -e '
  pins <- c(MESuSiE = "borangao/MESuSiE@5a3dd32",
            MACHINE = "xiangli428/MACHINE@f5abc5b",
            CARMA   = "ZikunY/CARMA@856f0fb")
  for (p in names(pins)) if (!requireNamespace(p, quietly = TRUE))
    remotes::install_github(pins[[p]], upgrade = "never", dependencies = FALSE)
  for (p in c("susieR", names(pins)))
    cat(sprintf("%-8s %s %s\n", p, packageVersion(p),
                if (is.null(packageDescription(p)$RemoteSha)) "(CRAN)" else packageDescription(p)$RemoteSha))'

# 3. SuSiEx (C++, built from source) and RSparsePro (a script), at pinned commits
mkdir -p "$TOOLS"
[ -d "$TOOLS/SuSiEx" ] || git clone -q https://github.com/getian107/SuSiEx "$TOOLS/SuSiEx"
git -C "$TOOLS/SuSiEx" checkout -q 1db2f58
[ -x "$TOOLS/SuSiEx/bin/SuSiEx" ] || (cd "$TOOLS/SuSiEx/src" && PATH="$PREFIX/bin:$PATH" make all)
[ -d "$TOOLS/RSparsePro_LD" ] || git clone -q https://github.com/zhwm/RSparsePro_LD "$TOOLS/RSparsePro_LD"
git -C "$TOOLS/RSparsePro_LD" checkout -q 5990896

# 4. FINEMAP v1.4.2, the official build every run used, md5-checked
FM="$TOOLS/finemap_v1.4.2/finemap_v1.4.2_x86_64"
if [ ! -x "$FM" ]; then
  mkdir -p "$TOOLS/finemap_v1.4.2"
  curl -sL -o /tmp/finemap_v1.4.2.tgz http://www.christianbenner.com/finemap_v1.4.2_x86_64.tgz
  tar xzf /tmp/finemap_v1.4.2.tgz -C "$TOOLS/finemap_v1.4.2" --strip-components=1
fi
echo "b5741e4374396696f511f2126e72efba  $FM" | md5sum -c -

# 5. record the fully resolved env next to the spec
mamba env export -p "$PREFIX" --no-builds | grep -v "^prefix:" > "$HERE/finemap_env.lock.yml"
echo "env ready: $PREFIX   tools: $TOOLS"

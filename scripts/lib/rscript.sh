#!/usr/bin/env bash
# Run Rscript from the environment recorded in project.conf.
#
# A conda R on Windows needs the environment's Library/bin on PATH for
# its runtime DLLs, which conda activation normally provides. Stages
# call this wrapper instead of a bare Rscript so that they work both
# with the conda environment built by 01_install.sh and with a system R.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
# shellcheck source=/dev/null
source "$ROOT/project.conf"
if [[ -n "${CONDA_PREFIX_PATH:-}" && -d "$CONDA_PREFIX_PATH" ]]; then
    export PATH="$CONDA_PREFIX_PATH/bin:$CONDA_PREFIX_PATH/Library/bin:$CONDA_PREFIX_PATH/Library/mingw-w64/bin:$CONDA_PREFIX_PATH/Scripts:$CONDA_PREFIX_PATH:$PATH"
    export CONDA_PREFIX="$CONDA_PREFIX_PATH"
    for r in "$CONDA_PREFIX_PATH/bin/Rscript" "$CONDA_PREFIX_PATH/Scripts/Rscript.exe" "$CONDA_PREFIX_PATH/Lib/R/bin/Rscript.exe"; do
        if [[ -x "$r" ]]; then exec "$r" "$@"; fi
    done
fi
exec "${RSCRIPT_SYSTEM:-Rscript}" "$@"

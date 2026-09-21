#!/usr/bin/env bash
# Build the software environment and record what was actually installed.
#
# Two routes, tried in order:
#   1. conda: create the msqrob2_tmt environment from config/environment.yml
#      (R >= 4.6, shellcheck, quarto), then let BiocManager install the
#      Bioconductor stack inside that R.
#   2. system R: if conda is absent or the environment fails to build,
#      use the Rscript already on PATH and install missing packages there.
#
# Either way the versions that ended up in use are written to
# config/r_packages.tsv and logs/sessionInfo_install.txt, and project.conf
# is updated to point at the Rscript that later stages must call.
#
# Usage: scripts/01_install.sh [--system-r] [--force] [--help]

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=/dev/null
source "$ROOT/project.conf"
LOG="$ROOT/logs/01_install.log"
STAMP="$ROOT/logs/01_install.done"
SYSTEM_R=0
FORCE=0

for arg in "$@"; do
    case "$arg" in
        --system-r) SYSTEM_R=1 ;;
        --force) FORCE=1 ;;
        --help|-h) sed -n '2,15p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "unknown argument: $arg" >&2; exit 2 ;;
    esac
done

mkdir -p "$ROOT/logs" "$ROOT/config"

if [[ -f "$STAMP" && $FORCE -eq 0 ]]; then
    echo "01_install: already completed on $(cat "$STAMP"); skipping (use --force to redo)."
    exit 0
fi

START=$(date +%s)

find_conda() {
    if command -v conda >/dev/null 2>&1; then command -v conda; return; fi
    for c in "$HOME/miniconda3/Scripts/conda.exe" "$HOME/miniconda3/bin/conda" \
             "$HOME/anaconda3/Scripts/conda.exe" "$HOME/anaconda3/bin/conda" \
             "$HOME/miniforge3/bin/conda" "$HOME/mambaforge/bin/conda" \
             "/opt/conda/bin/conda"; do
        if [[ -x "$c" ]]; then echo "$c"; return; fi
    done
    echo ""
}

set_conf() {
    # Replace KEY=... in project.conf (portable sed -i, no backup file).
    local key="$1" val="$2"
    if grep -q "^$key=" "$ROOT/project.conf"; then
        sed -i.bak "s|^$key=.*|$key=\"$val\"|" "$ROOT/project.conf" && rm -f "$ROOT/project.conf.bak"
    else
        echo "$key=\"$val\"" >> "$ROOT/project.conf"
    fi
}

CONDA_BIN=""
if [[ $SYSTEM_R -eq 0 ]]; then CONDA_BIN="$(find_conda)"; fi

ENV_FOUND=0
if [[ -n "$CONDA_BIN" ]]; then
    echo "01_install: conda found at $CONDA_BIN"
    ENV_NAME="msqrob2_tmt"
    if "$CONDA_BIN" env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
        echo "01_install: conda environment $ENV_NAME already exists"
    else
        echo "01_install: creating conda environment $ENV_NAME (this downloads R and takes minutes)"
        if ! "$CONDA_BIN" env create -f "$ROOT/config/environment.yml" -n "$ENV_NAME" -y 2>&1 | tee -a "$LOG"; then
            echo "01_install: conda environment creation failed; falling back to system R" | tee -a "$LOG"
            CONDA_BIN=""
        fi
    fi
    if [[ -n "$CONDA_BIN" ]]; then
        PREFIX="$("$CONDA_BIN" env list | awk -v n="$ENV_NAME" '$1==n {print $NF}' | tr -d '\r')"
        # Git Bash needs a POSIX path; conda prints a Windows path there.
        if command -v cygpath >/dev/null 2>&1; then PREFIX="$(cygpath -u "$PREFIX")"; fi
        if [[ -x "$PREFIX/bin/Rscript" || -x "$PREFIX/Scripts/Rscript.exe" ]]; then
            ENV_FOUND=1
            set_conf CONDA_ENV "$ENV_NAME"
            set_conf CONDA_EXE_PATH "$CONDA_BIN"
            set_conf CONDA_PREFIX_PATH "$PREFIX"
            # Keep the conda quarto only if it starts; the conda-forge
            # launcher on Windows was broken when this was written, and
            # 00_configure.sh may already have found a working one.
            for q in "$PREFIX/bin/quarto" "$PREFIX/Library/bin/quarto.cmd" "$PREFIX/Scripts/quarto.cmd"; do
                if [[ -e "$q" ]] && "$q" --version >/dev/null 2>&1; then set_conf QUARTO "$q"; break; fi
            done
            for s in "$PREFIX/bin/shellcheck" "$PREFIX/Library/bin/shellcheck.exe" "$PREFIX/Scripts/shellcheck.exe"; do
                if [[ -x "$s" ]]; then set_conf SHELLCHECK "$s"; break; fi
            done
        else
            echo "01_install: could not locate Rscript inside $PREFIX; falling back to system R" | tee -a "$LOG"
        fi
    fi
fi

# Stages never call Rscript directly. scripts/lib/rscript.sh reads
# project.conf and either puts the conda environment's runtime
# directories on PATH (a conda R on Windows cannot start without them)
# or falls through to the system Rscript recorded in RSCRIPT_SYSTEM.
set_conf RSCRIPT_SYSTEM "${RSCRIPT:-$(command -v Rscript || true)}"
if [[ $ENV_FOUND -eq 0 ]]; then
    set_conf CONDA_PREFIX_PATH ""
    if [[ -z "$(command -v Rscript || true)" ]]; then
        echo "01_install: no Rscript found and conda unavailable; install R >= 4.4 and re-run" >&2
        exit 1
    fi
    echo "01_install: using system R at $(command -v Rscript)"
fi
RSCRIPT_USE="$ROOT/scripts/lib/rscript.sh"

# Install the Bioconductor stack inside whichever R was chosen. On a
# conda R for Windows, BiocManager fetches CRAN and Bioconductor
# binaries built for the same UCRT toolchain; on Linux it compiles.
echo "01_install: installing R packages with $("$RSCRIPT_USE" --version 2>&1 | head -1)" | tee -a "$LOG"
if ! "$RSCRIPT_USE" "$ROOT/scripts/lib/install_packages.R" \
        "$ROOT/config/r_packages.tsv" "$ROOT/logs/sessionInfo_install.txt" 2>&1 | tee -a "$LOG"; then
    if [[ $ENV_FOUND -eq 1 ]]; then
        echo "01_install: package installation failed in the conda R; retrying with system R" | tee -a "$LOG"
        set_conf CONDA_ENV ""
        set_conf CONDA_PREFIX_PATH ""
        "$RSCRIPT_USE" "$ROOT/scripts/lib/install_packages.R" \
            "$ROOT/config/r_packages.tsv" "$ROOT/logs/sessionInfo_install.txt" 2>&1 | tee -a "$LOG"
    else
        exit 1
    fi
fi

set_conf RSCRIPT "$RSCRIPT_USE"

# shellcheck is optional for running the pipeline but the final checks use it.
if ! grep -q '^SHELLCHECK=' "$ROOT/project.conf"; then
    set_conf SHELLCHECK "$(command -v shellcheck || true)"
fi

END=$(date +%s)
printf "stage\telapsed_s\tpeak_rss_mb\tnote\n01_install\t%d\tNA\tRscript=%s\n" "$((END-START))" "$RSCRIPT_USE" > "$ROOT/logs/01_install.resources.tsv"
date -u +%Y-%m-%dT%H:%M:%SZ > "$STAMP"
echo "01_install: done in $((END-START)) s; R packages recorded in config/r_packages.tsv"

#!/usr/bin/env bash
# Detect or accept the machine's CPU, RAM and disk budget and write
# project.conf at the repository root. Every later stage sources that
# file instead of hardcoding resources, so the same scripts run on a
# 4 GB laptop and on a workstation without edits.
#
# Usage:
#   scripts/00_configure.sh [--threads N] [--ram GB] [--disk GB] [--yes] [--force] [--help]
#
# Without arguments the script detects the values and asks for
# confirmation; --yes accepts the detected values non-interactively.
# Re-running with an existing project.conf and no --force skips.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONF="$ROOT/project.conf"
THREADS=""
RAM_GB=""
DISK_GB=""
YES=0
FORCE=0
RAM_GB_EXPLICIT=""

usage() {
    sed -n '2,12p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --threads) THREADS="$2"; shift 2 ;;
        --ram)     RAM_GB="$2"; RAM_GB_EXPLICIT=1; shift 2 ;;
        --disk)    DISK_GB="$2"; shift 2 ;;
        --yes)     YES=1; shift ;;
        --force)   FORCE=1; shift ;;
        --help|-h) usage; exit 0 ;;
        *) echo "unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

if [[ -f "$CONF" && $FORCE -eq 0 && -z "$THREADS$RAM_GB$DISK_GB" ]]; then
    echo "project.conf already exists; skipping configuration (use --force to rewrite)."
    exit 0
fi

# Detection. The three platforms this has to work on are Linux, macOS
# and Git Bash on Windows; the last one has no /proc/meminfo and its
# nproc reports logical processors, which is what we want anyway.
detect_threads() {
    if command -v nproc >/dev/null 2>&1; then nproc
    elif command -v sysctl >/dev/null 2>&1; then sysctl -n hw.ncpu
    else echo 2; fi
}

detect_ram_gb() {
    if [[ -r /proc/meminfo ]]; then
        awk '/MemTotal/ {printf "%d", $2/1024/1024}' /proc/meminfo
    elif command -v sysctl >/dev/null 2>&1 && sysctl -n hw.memsize >/dev/null 2>&1; then
        echo $(( $(sysctl -n hw.memsize) / 1024 / 1024 / 1024 ))
    elif command -v powershell.exe >/dev/null 2>&1; then
        powershell.exe -NoProfile -Command \
            "[int][math]::Floor((Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory/1GB)" \
            | tr -d '\r'
    else
        echo 4
    fi
}

detect_disk_gb() {
    # Free space on the filesystem holding the repository, in whole GB.
    df -Pk "$ROOT" | awk 'NR==2 {printf "%d", $4/1024/1024}'
}

DET_THREADS="$(detect_threads)"
DET_RAM="$(detect_ram_gb)"
DET_DISK="$(detect_disk_gb)"

THREADS="${THREADS:-$DET_THREADS}"
RAM_GB="${RAM_GB:-$DET_RAM}"
DISK_GB="${DISK_GB:-$DET_DISK}"

# The RAM figure that the R stages plan against. The brief sets a hard
# ceiling of roughly 4 GB usable, so we plan for min(detected, 4) unless
# the user passes --ram explicitly. Detected total RAM is not usable RAM:
# on the build machine 16 GB were installed and 1.3 GB were free at the
# time of detection because of other processes.
PLAN_RAM_GB="$RAM_GB"
if [[ -z "$RAM_GB_EXPLICIT" ]]; then
    if [[ "$RAM_GB" -gt 4 ]]; then PLAN_RAM_GB=4; fi
fi

# BiocParallel workers. Serial by default: the PSM-level mixed models
# copy the QFeatures object into each worker, and at 4 GB two copies of
# a quarter-gigabyte PSM table plus lme4 working memory exhaust RAM
# before they save wall clock. run_all.sh --workers N overrides this.
WORKERS=1

# Locate R and quarto. Later scripts call "$RSCRIPT" and "$QUARTO" so
# that a conda environment created by 01_install.sh can be preferred
# over a system R without editing every script.
R_CANDIDATE="$(command -v Rscript || true)"
QUARTO_CANDIDATE="$(command -v quarto || true)"
if [[ -z "$QUARTO_CANDIDATE" ]]; then
    for q in "/c/Program Files/RStudio/resources/app/bin/quarto/bin/quarto.exe" \
             "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto" \
             "/usr/lib/rstudio/resources/app/bin/quarto/bin/quarto"; do
        if [[ -x "$q" ]]; then QUARTO_CANDIDATE="$q"; break; fi
    done
fi

cat <<EOF
Detected resources
  threads : $DET_THREADS (using $THREADS)
  RAM     : $DET_RAM GB installed (planning against $PLAN_RAM_GB GB)
  disk    : $DET_DISK GB free (using $DISK_GB)
  Rscript : ${R_CANDIDATE:-not found}
  quarto  : ${QUARTO_CANDIDATE:-not found}
EOF

if [[ $YES -eq 0 ]]; then
    read -r -p "Write these values to project.conf? [y/N] " ans
    case "$ans" in y|Y|yes|YES) ;; *) echo "aborted"; exit 1 ;; esac
fi

cat > "$CONF" <<EOF
# Written by scripts/00_configure.sh on $(date -u +%Y-%m-%dT%H:%M:%SZ)
# Source this file; do not edit by hand, re-run the script instead.
PROJECT_ROOT="$ROOT"
THREADS=$THREADS
RAM_GB=$RAM_GB
PLAN_RAM_GB=$PLAN_RAM_GB
DISK_GB=$DISK_GB
WORKERS=$WORKERS
RSCRIPT="${R_CANDIDATE}"
QUARTO="${QUARTO_CANDIDATE}"
# Set by 01_install.sh when a conda environment provides R
CONDA_ENV=""
CONDA_EXE_PATH=""
EOF

echo "wrote $CONF"

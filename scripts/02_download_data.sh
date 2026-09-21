#!/usr/bin/env bash
# Fetch and verify the input tables.
#
# Zenodo record 14767905 is enumerated through its API; files are chosen
# by dataset prefix, downloaded once into a BiocFileCache under
# data/.bfc, verified against the record's md5, and linked into data/.
# The MS2-only counterpart of the spike-in is fetched from MassIVE
# MSV000084266 the same way, with its size checked against the MassIVE
# listing and its md5 recorded at first retrieval.
#
# Usage: scripts/02_download_data.sh [--dataset spikein2] [--no-ms2] [--force] [--help]
#
# Default datasets: spikein1 and mouse (about 320 MB from Zenodo) plus the
# MS2 PSM export (288 MB from MassIVE). --dataset spikein2 adds the
# msTrawler benchmark tables (24 MB).

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=/dev/null
source "$ROOT/project.conf"
STAMP="$ROOT/logs/02_download_data.done"
DATASETS="spikein1,mouse"
EXTRA_ARGS=()
FORCE=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dataset)
            if [[ "$2" == "spikein2" ]]; then DATASETS="$DATASETS,spikein2"; else echo "unknown dataset $2" >&2; exit 2; fi
            shift 2 ;;
        --no-ms2) EXTRA_ARGS+=("--no-ms2"); shift ;;
        --force) FORCE=1; EXTRA_ARGS+=("--force"); shift ;;
        --help|-h) sed -n '2,15p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "unknown argument: $1" >&2; exit 2 ;;
    esac
done

mkdir -p "$ROOT/logs" "$ROOT/data"

# The stamp records which datasets were fetched, so asking for spikein2
# later re-runs the stage instead of skipping it.
if [[ -f "$STAMP" && $FORCE -eq 0 ]] && grep -qx "datasets=$DATASETS" "$STAMP"; then
    echo "02_download_data: already completed for $DATASETS on $(head -1 "$STAMP"); skipping (use --force to redo)."
    exit 0
fi

# Disk check before touching the network: the default fetch needs about
# 0.65 GB, spikein2 adds 0.03 GB, and the cache plus links share inodes.
FREE_GB=$(df -Pk "$ROOT" | awk 'NR==2 {printf "%d", $4/1024/1024}')
if [[ "$FREE_GB" -lt 2 ]]; then
    echo "02_download_data: only ${FREE_GB} GB free; need at least 2 GB" >&2
    exit 1
fi

START=$(date +%s)
"$RSCRIPT" "$ROOT/scripts/lib/download_data.R" --datasets "$DATASETS" "${EXTRA_ARGS[@]+"${EXTRA_ARGS[@]}"}"
END=$(date +%s)

printf "stage\telapsed_s\tpeak_rss_mb\tnote\n02_download_data\t%d\tNA\tdatasets=%s\n" \
    "$((END-START))" "$DATASETS" > "$ROOT/logs/02_download_data.resources.tsv"
{ date -u +%Y-%m-%dT%H:%M:%SZ; echo "datasets=$DATASETS"; } > "$STAMP"
echo "02_download_data: done in $((END-START)) s"

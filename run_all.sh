#!/usr/bin/env bash
# Run the whole analysis in order. Every stage is idempotent and skips
# itself when its stamp file in logs/ exists, so the script can be
# re-run after a failure and picks up where it stopped.
#
# Usage:
#   bash run_all.sh [--from <stage>] [--dataset spikein2] [--search-arm]
#                   [--workers N] [--psm-budget-min M] [--force] [--help]
#
#   --from <stage>      start at a stage number (00 to 10) or name prefix,
#                       for example --from 05
#   --dataset spikein2  also fetch the msTrawler multibatch tables (not
#                       analysed by the core stages; see README)
#   --search-arm        run the optional Sage re-search of one mixture
#                       (scripts/10_search_arm.sh); disk-gated
#   --workers N         BiocParallel workers for the model stages; the
#                       default is 1 because workers copy the data
#   --psm-budget-min M  wall clock budget for the PSM-level spike-in fits
#   --force             re-run stages even if their stamp exists
#
# Resource use per stage is written to logs/<stage>.resources.tsv.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FROM="00"
DATASET_ARGS=()
SEARCH_ARM=0
WORKERS=""
PSM_BUDGET=""
FORCE=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --from) FROM="$2"; shift 2 ;;
        --dataset) DATASET_ARGS+=(--dataset "$2"); shift 2 ;;
        --search-arm) SEARCH_ARM=1; shift ;;
        --workers) WORKERS="$2"; shift 2 ;;
        --psm-budget-min) PSM_BUDGET="$2"; shift 2 ;;
        --force) FORCE=(--force); shift ;;
        --help|-h) sed -n '2,21p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "unknown argument: $1" >&2; exit 2 ;;
    esac
done

FROM_NUM="${FROM:0:2}"
runs_from() { [[ "$1" -ge "$FROM_NUM" ]]; }

mkdir -p "$ROOT/logs" "$ROOT/results" "$ROOT/figures" "$ROOT/data"
T_ALL=$(date +%s)

stage() {
    # stage <number> <description> <command...>
    local num="$1" desc="$2"; shift 2
    if ! runs_from "$num"; then return 0; fi
    echo "==> [$num] $desc"
    local t0; t0=$(date +%s)
    "$@" 2>&1 | tee "$ROOT/logs/${num}_console.log"
    local rc=${PIPESTATUS[0]}
    if [[ $rc -ne 0 ]]; then echo "stage $num failed with exit code $rc" >&2; exit "$rc"; fi
    echo "    [$num] done in $(( $(date +%s) - t0 )) s"
}

stage 00 "configure" bash "$ROOT/scripts/00_configure.sh" --yes
# shellcheck source=/dev/null
source "$ROOT/project.conf"
stage 01 "install" bash "$ROOT/scripts/01_install.sh" "${FORCE[@]+"${FORCE[@]}"}"
# project.conf may now point at a conda environment
# shellcheck source=/dev/null
source "$ROOT/project.conf"
R="$ROOT/scripts/lib/rscript.sh"

WORKER_ARGS=()
if [[ -n "$WORKERS" ]]; then WORKER_ARGS=(--workers "$WORKERS"); fi
BUDGET_ARGS=()
if [[ -n "$PSM_BUDGET" ]]; then BUDGET_ARGS=(--psm-budget-min "$PSM_BUDGET"); fi

stage 02 "download and verify data" bash "$ROOT/scripts/02_download_data.sh" "${DATASET_ARGS[@]+"${DATASET_ARGS[@]}"}" "${FORCE[@]+"${FORCE[@]}"}"
stage 03 "build QFeatures objects" "$R" "$ROOT/scripts/03_build_qfeatures.R" "${FORCE[@]+"${FORCE[@]}"}"
stage 04 "preprocess" "$R" "$ROOT/scripts/04_preprocess.R" "${FORCE[@]+"${FORCE[@]}"}"
stage 05 "spike-in benchmark" "$R" "$ROOT/scripts/05_benchmark_workflows.R" "${WORKER_ARGS[@]+"${WORKER_ARGS[@]}"}" "${BUDGET_ARGS[@]+"${BUDGET_ARGS[@]}"}" "${FORCE[@]+"${FORCE[@]}"}"
stage 06 "ratio compression" "$R" "$ROOT/scripts/06_compression.R" "${FORCE[@]+"${FORCE[@]}"}"
stage 07 "mouse inference" "$R" "$ROOT/scripts/07_mouse_inference.R" "${WORKER_ARGS[@]+"${WORKER_ARGS[@]}"}" "${FORCE[@]+"${FORCE[@]}"}"
stage 08 "figures" "$R" "$ROOT/scripts/08_figures.R" "${FORCE[@]+"${FORCE[@]}"}"
stage 09 "report" bash "$ROOT/scripts/09_render_report.sh" "${FORCE[@]+"${FORCE[@]}"}"

if [[ $SEARCH_ARM -eq 1 ]]; then
    stage 10 "Sage search arm" bash "$ROOT/scripts/10_search_arm.sh" "${FORCE[@]+"${FORCE[@]}"}"
else
    echo "==> [10] Sage search arm not requested (pass --search-arm); core results do not depend on it"
fi

# One table with every stage's measured cost, for the README.
{
    printf "stage\telapsed_s\tpeak_rss_mb\tnote\n"
    for f in "$ROOT"/logs/*.resources.tsv; do
        [[ -f "$f" ]] && tail -n +2 "$f" | cut -f1-4
    done
} > "$ROOT/logs/resources_all.tsv"
echo "==> all stages finished in $(( $(date +%s) - T_ALL )) s; per-stage cost in logs/resources_all.tsv"

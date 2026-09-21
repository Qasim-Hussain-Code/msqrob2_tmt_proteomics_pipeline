#!/usr/bin/env bash
# Render scripts/09_report.qmd to results/report.html with the R that
# the pipeline uses. Quarto picks its R from QUARTO_R, so the wrapper
# resolves the same binary that scripts/lib/rscript.sh runs.
#
# Usage: scripts/09_render_report.sh [--force]

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=/dev/null
source "$ROOT/project.conf"
STAMP="$ROOT/logs/09_report.done"
FORCE=0
for arg in "$@"; do
    case "$arg" in
        --force) FORCE=1 ;;
        --help|-h) sed -n '2,6p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
    esac
done
if [[ -f "$STAMP" && $FORCE -eq 0 ]]; then
    echo "09_report: already completed on $(cat "$STAMP"); skipping (use --force to redo)."
    exit 0
fi
if [[ -z "${QUARTO:-}" || ! -e "$QUARTO" ]]; then
    echo "09_report: quarto not found (QUARTO in project.conf); install quarto or set QUARTO" >&2
    exit 3
fi

# The stage resource table the report reads; run_all.sh rebuilds it at
# the end as well, so this one may miss later stages.
{
    printf "stage\telapsed_s\tpeak_rss_mb\tnote\n"
    for f in "$ROOT"/logs/*.resources.tsv; do
        [[ -f "$f" ]] && tail -n +2 "$f" | cut -f1-4
    done
} > "$ROOT/logs/resources_all.tsv"

START=$(date +%s)
if [[ -n "${CONDA_PREFIX_PATH:-}" && -d "$CONDA_PREFIX_PATH" ]]; then
    export PATH="$CONDA_PREFIX_PATH/bin:$CONDA_PREFIX_PATH/Library/bin:$CONDA_PREFIX_PATH/Library/mingw-w64/bin:$CONDA_PREFIX_PATH/Scripts:$PATH"
    for r in "$CONDA_PREFIX_PATH/bin/R" "$CONDA_PREFIX_PATH/Scripts/R.exe" "$CONDA_PREFIX_PATH/Lib/R/bin/R.exe"; do
        if [[ -x "$r" ]]; then export QUARTO_R="$r"; break; fi
    done
fi
mkdir -p "$ROOT/results"
( cd "$ROOT/scripts" && "$QUARTO" render 09_report.qmd --to html --output report.html --output-dir "$ROOT/results" )
# Quarto embeds the whole Bootstrap Icons stylesheet; one icon rule is
# removed afterwards, see scripts/lib/strip_icon_rule.R.
"$RSCRIPT" "$ROOT/scripts/lib/strip_icon_rule.R" "$ROOT/results/report.html"
END=$(date +%s)
printf "stage\telapsed_s\tpeak_rss_mb\tnote\n09_report\t%d\tNA\tquarto %s\n" "$((END-START))" "$("$QUARTO" --version 2>/dev/null | head -1)" > "$ROOT/logs/09_report.resources.tsv"
date -u +%Y-%m-%dT%H:%M:%SZ > "$STAMP"
echo "09_report: rendered results/report.html in $((END-START)) s"

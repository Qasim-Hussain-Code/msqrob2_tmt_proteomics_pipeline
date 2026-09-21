#!/usr/bin/env bash
# The checks run before the repository is handed over. Each prints what
# it found; the exit status is 1 if any check fails.
#
# Usage: scripts/lib/final_checks.sh
set -uo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT" || exit 1
# shellcheck source=/dev/null
source "$ROOT/project.conf"
status=0

echo "== 1. shellcheck on every shell script =="
if [[ -n "${SHELLCHECK:-}" && -x "$SHELLCHECK" ]]; then
    if "$SHELLCHECK" -x scripts/*.sh scripts/lib/*.sh run_all.sh; then echo "clean"; else status=1; fi
else
    echo "shellcheck not available"; status=1
fi

echo "== 2. every figure the README embeds exists and is written by a script =="
while read -r fig; do
    f="figures/$fig"
    if [[ ! -s "$f" ]]; then echo "MISSING: $f"; status=1; continue; fi
    if ! grep -q -- "$fig" scripts/*.R; then
        # figures named by pattern (dataset prefix) are matched on their suffix
        suffix="${fig#*_}"
        if ! grep -q -- "$suffix" scripts/*.R; then echo "NOT PRODUCED BY A SCRIPT: $f"; status=1; fi
    fi
done < <(grep -oE 'figures/[A-Za-z0-9_./-]+\.png' README.md | sed 's|figures/||' | sort -u)
echo "checked $(grep -oE 'figures/[A-Za-z0-9_./-]+\.png' README.md | sort -u | wc -l) figures"

echo "== 3. no em dash, emoji or banned phrase =="
if bash scripts/lib/check_prose.sh > /tmp/prose_check.txt; then echo "clean"; else cat /tmp/prose_check.txt; status=1; fi

echo "== 4. no tracked file over 50 MB =="
big=$(git ls-files -z | xargs -0 stat -c '%s %n' 2>/dev/null | awk '$1 > 52428800' | sort -n)
if [[ -n "$big" ]]; then echo "$big"; status=1; else echo "largest tracked file: $(git ls-files -z | xargs -0 stat -c '%s %n' | sort -n | tail -1)"; fi

echo "== 5. model formulas match in the three stated places =="
check_formula() {
    local label="$1" formula="$2" script="$3"
    local n_script n_readme
    n_script=$(grep -cF -- "$formula" "$script")
    n_readme=$(grep -cF -- "$formula" README.md)
    if [[ $n_script -ge 2 && $n_readme -ge 1 ]]; then
        echo "ok  $label: $n_script in $(basename "$script"), $n_readme in README"
    else
        echo "BAD $label: $n_script in $(basename "$script"), $n_readme in README"; status=1
    fi
}
check_formula "spike-in mixed"   '~ 0 + Condition + (1 | Run) + (1 | Mixture)' scripts/05_benchmark_workflows.R
check_formula "spike-in PSM"     '~ 0 + Condition + (1 | Run) + (1 | Mixture) + (1 | Run:Channel) + (1 | Run:ionID)' scripts/05_benchmark_workflows.R
check_formula "mouse fraction"   '~ Diet * Duration + (1 | Mixture) + (1 | Run) + (1 | BioReplicate)' scripts/07_mouse_inference.R
check_formula "mouse mixture"    '~ Diet * Duration + (1 | Mixture)' scripts/07_mouse_inference.R
check_formula "mouse PSM"        '~ Diet * Duration + (1 | Mixture) + (1 | Run) + (1 | Run:Channel) + (1 | Run:ionID) + (1 | BioReplicate)' scripts/07_mouse_inference.R

echo "== 6. stage stamps and resource logs =="
for st in 01_install 02_download_data 03_build_qfeatures 04_preprocess 05_benchmark_workflows 06_compression 07_mouse_inference 08_figures 09_report; do
    if [[ -f "logs/$st.resources.tsv" ]]; then printf "%-24s %s\n" "$st" "$(tail -1 "logs/$st.resources.tsv" | cut -f2,3)"; else echo "$st: no resource log"; status=1; fi
done

exit $status

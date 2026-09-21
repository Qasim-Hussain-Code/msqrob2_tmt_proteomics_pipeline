#!/usr/bin/env bash
# Search the repository's prose and code for em dashes, emoji and the
# phrases the writing brief bans. Exit 1 if anything is found.
#
# Usage: scripts/lib/check_prose.sh
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT"

# Tracked text files only: data/ and results/qfeatures are never prose.
mapfile -t FILES < <(git ls-files | grep -vE '\.(png|gz|rds|pdf|ico|html)$')

status=0
echo "== em dashes (U+2014) =="
if grep -nP '\x{2014}' "${FILES[@]}" 2>/dev/null; then status=1; else echo "none"; fi

echo "== emoji and pictographs (U+1F000 and above, U+2600 to U+27BF) =="
if grep -nP '[\x{1F000}-\x{1FFFF}\x{2600}-\x{27BF}]' "${FILES[@]}" 2>/dev/null; then status=1; else echo "none"; fi

echo "== banned phrases =="
PHRASES=(
    "it is worth noting" "it is important to note" "rapidly evolving" "plays a crucial role"
    "serves as a testament" "paving the way" "in conclusion" "delve" "leverag" "seamless"
    "comprehensive" "underscore" "showcase" "highlight" "not only" "it is not .* it is"
    "overall," "^overall " "robust framework" "robust approach" "robust pipeline"
)
for ph in "${PHRASES[@]}"; do
    if grep -niE -- "$ph" "${FILES[@]}" 2>/dev/null | grep -vE "check_prose.sh|robustSummary|robust regression|robust = |robust=|robust M-|M-estimat|\"robust\"|robust_|_robust|robust ridge|robust linear|robust fit|robust lmer|Robust Linear|Robust Summar|robust Peptide|Robust Ridge|robust weights|not robust|robust\)" ; then status=1; fi
done
if [[ $status -eq 0 ]]; then echo "no hits"; fi
exit $status

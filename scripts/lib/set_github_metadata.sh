#!/usr/bin/env bash
# Set the GitHub description and topics once the repository has a remote.
# Not run by the pipeline; run it yourself after `git push`.
#
# Usage: scripts/lib/set_github_metadata.sh [owner/repo]
set -euo pipefail
REPO="${1:-}"
ARGS=()
if [[ -n "$REPO" ]]; then ARGS=(--repo "$REPO"); fi
DESCRIPTION="TMT proteomics differential abundance with msqrob2 mixed models, benchmarked on spike-in mixtures with known ground truth to measure ratio compression and the cost of ignoring the plex, and applied to a high fat diet study in mouse adipose tissue."
TOPICS=(proteomics mass-spectrometry tmt isobaric-labeling tandem-mass-tags msqrob2 qfeatures bioconductor
        mixed-models differential-abundance ratio-compression reporter-ions sps-ms3 reference-channels
        spike-in-benchmark multiplexing batch-effects proteome-discoverer bioinformatics r)
gh repo edit "${ARGS[@]+"${ARGS[@]}"}" --description "$DESCRIPTION"
for t in "${TOPICS[@]}"; do
    gh repo edit "${ARGS[@]+"${ARGS[@]}"}" --add-topic "$t"
done
echo "description and ${#TOPICS[@]} topics set"

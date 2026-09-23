#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Every selectable method, run on the demo as far as superclusters.
#
#   bash tests/test_config_paths.sh
#
# Exists because auto_label: shape_simple was broken for its entire life and
# nothing noticed: it is not on the default path and the demo never selects
# it. A config option that nothing runs is an option that does not work.
#
# Separate from tier 2 because it is seven pipeline runs. It runs nightly in
# CI next to tier 3. Each run caps max_features hard: these assert that a
# path runs and produces a sane table, not what the table says.
# ---------------------------------------------------------------------------

source "$(dirname "${BASH_SOURCE[0]}")/lib/assert.sh"
cd "$REPO_ROOT" || exit 1

if [[ -n "${CONDA_PREFIX:-}" && -x "${CONDA_PREFIX}/bin/Rscript" ]]; then
    RUN=(); else RUN=(pixi run --frozen)
fi

SANDBOX="$(mktemp -d "${TMPDIR:-/tmp}/timecourse-patterns-paths-XXXXXX")"
cleanup() {
    if [[ "${KEEP_SANDBOX:-0}" == "1" ]] || (( FAIL > 0 )); then
        printf "\nSandbox kept: %s\n" "$SANDBOX"
    else
        rm -rf "$SANDBOX"
    fi
}
trap cleanup EXIT

# name | override YAML. Every override is deep-merged onto config/demo.yaml.
run_path () {
    local name="$1" override="$2" dir="$SANDBOX/$1"
    mkdir -p "$dir"
    printf "%s\noutput:\n  dir: %s/results\n" "$override" "$dir" > "$dir/override.yaml"
    "${RUN[@]}" python tests/helpers/merge_config.py config/demo.yaml "$dir/override.yaml" \
        > "$dir/config.yaml"
    if "${RUN[@]}" snakemake --configfile "$dir/config.yaml" -j 4 \
            --until superclusters > "$dir/run.log" 2>&1; then
        local t n
        t=$(find "$dir/results/clusters" -name '*_clusters.tsv' | head -1)
        n=$( [[ -n "$t" ]] && tail -n +2 "$t" | wc -l || echo 0)
        if (( n > 0 )); then pass "$name ($n features assigned)"
        else fail "$name completed but its cluster table is empty"; fi
    else
        fail "$name failed (see $dir/run.log)"
        grep -E 'Error|Exception' "$dir/run.log" | head -3
    fi
}

group "Alternative config paths"
CAP="cluster:
  max_features: 600"

run_path wald_union   "differential:
  method: wald_union
  min_contrasts: 2
$CAP"
run_path method_none  "differential:
  method: none
  top_n: 2000
$CAP"
run_path shape_simple "superclusters:
  auto_label: shape_simple
  split: null
$CAP"
run_path hclust       "superclusters:
  method: hclust
  k: 3
  split: null
$CAP"
run_path minc_auto    "cluster:
  max_features: 600
  degpatterns:
    minc: auto"
run_path subsample    "cluster:
  max_features: 400"
run_path gene_mode    "input:
  features: null
$CAP"

# Gene mode must not write coordinates or BED files.
g="$SANDBOX/gene_mode/results/clusters/WT_clusters.tsv"
if [[ -f "$g" ]]; then
    if head -1 "$g" | grep -q 'chr'; then fail "gene mode wrote coordinate columns"
    else pass "gene mode cluster table has no coordinate columns"; fi
fi

summary

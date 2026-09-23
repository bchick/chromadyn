#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Tier 2: run the bundled demo end to end and check what came out.
#
#   bash tests/test_endtoend.sh
#
# Runs in a sandbox so it cannot disturb a working results/ directory. The
# sandbox is deleted on success and PRESERVED on failure, with its path
# printed, because the first thing you want after a failure is the output.
# Set KEEP_SANDBOX=1 to keep it either way.
#
# Golden tables are diffed against tests/golden/. When a change to the method
# is intended, regenerate them with UPDATE_GOLDEN=1.
# ---------------------------------------------------------------------------

source "$(dirname "${BASH_SOURCE[0]}")/lib/assert.sh"
cd "$REPO_ROOT" || exit 1

if [[ -n "${CONDA_PREFIX:-}" && -x "${CONDA_PREFIX}/bin/Rscript" ]]; then
    RUN=(); else RUN=(pixi run --frozen)
fi

SANDBOX="$(mktemp -d "${TMPDIR:-/tmp}/timecourse-patterns-e2e-XXXXXX")"
CFG="$SANDBOX/config.yaml"
OUT="$SANDBOX/results"
GOLDEN="$REPO_ROOT/tests/golden"

cleanup() {
    if [[ "${KEEP_SANDBOX:-0}" == "1" ]]; then
        printf "\n%sSandbox kept:%s %s\n" "$CYAN" "$RESET" "$SANDBOX"
    elif (( FAIL > 0 )); then
        printf "\n%sSandbox preserved for inspection:%s %s\n" "$YELLOW" "$RESET" "$SANDBOX"
    else
        rm -rf "$SANDBOX"
    fi
}
trap cleanup EXIT

sed "s#^  dir: results#  dir: $OUT#" config/demo.yaml > "$CFG"

# ---------------------------------------------------------------------------
group "Full run on the bundled demo"
# ---------------------------------------------------------------------------
if "${RUN[@]}" snakemake --configfile "$CFG" -j 4 > "$SANDBOX/run.log" 2>&1; then
    pass "pipeline completed"
else
    fail "pipeline failed (see $SANDBOX/run.log)"
    tail -30 "$SANDBOX/run.log"
    summary; exit 1
fi

ARMS=$(awk -F'\t' 'NR>1 {print $2}' demo/samplesheet.tsv | sort -u | grep -v '^shared$')

# ---------------------------------------------------------------------------
group "Every declared output exists and is non-empty"
# ---------------------------------------------------------------------------
for f in qc/validation_report.tsv qc/library_sizes.tsv qc/replicate_correlation.tsv \
         clusters/cross_arm_crosstab.tsv objects/dds.rds objects/transformed.rds \
         report.html run_manifest.json; do
    assert_file "$f" "$OUT/$f"
done
for fmt in pdf svg png; do
    for f in pca sample_correlation replicate_correlation; do
        assert_file "qc/$f.$fmt" "$OUT/qc/$f.$fmt"
    done
done
for arm in $ARMS; do
    for f in "differential/${arm}_results.tsv.gz" \
             "clusters/${arm}_clusters.tsv" "clusters/${arm}_cluster_profiles.tsv" \
             "clusters/${arm}_supercluster_sizes.tsv" "clusters/${arm}_k_diagnostics.tsv" \
             "clusters/${arm}_cluster_assignment_qc.tsv" \
             "objects/${arm}_differential.rds" "objects/${arm}_degpatterns.rds"; do
        assert_file "$f" "$OUT/$f"
    done
    for fmt in pdf svg png; do
        for f in dynamic_vs_static degpatterns dendrogram supercluster_traces \
                 supercluster_ribbon heatmap; do
            assert_file "figures/${arm}_${f}.${fmt}" "$OUT/figures/${arm}_${f}.${fmt}"
        done
    done
done

# ---------------------------------------------------------------------------
group "Cluster table schema and content"
# ---------------------------------------------------------------------------
for arm in $ARMS; do
    t="$OUT/clusters/${arm}_clusters.tsv"
    hdr=$(head -1 "$t")
    assert_eq "cluster table header ($arm)" "$hdr" \
        "feature_id	cluster	supercluster	supercluster_label	assign_method	assign_cor	chr	start	end"
    n=$(tail -n +2 "$t" | wc -l)
    if (( n > 0 )); then pass "cluster table has $n rows ($arm)"; else fail "cluster table empty ($arm)"; fi
    nsc=$(tail -n +2 "$t" | cut -f4 | sort -u | wc -l)
    if (( nsc >= 2 && nsc <= 12 )); then
        pass "distinct trajectory classes in range ($arm): $nsc"
    else
        fail "distinct trajectory classes out of range ($arm): $nsc, expected 2 to 12"
    fi
    # No feature may appear twice.
    dup=$(tail -n +2 "$t" | cut -f1 | sort | uniq -d | wc -l)
    assert_eq "no duplicated feature_id ($arm)" "$dup" "0"
done

# ---------------------------------------------------------------------------
group "BED exports (guardrail 9.4)"
# ---------------------------------------------------------------------------
# The source analysis silently wrote thousands of rows of "chrN.start.end NA NA".
shopt -s nullglob
beds=("$OUT"/clusters/bed/*/*.bed)
if (( ${#beds[@]} == 0 )); then
    fail "no BED files written in region mode"
else
    pass "${#beds[@]} BED files written"
    bad=0
    for b in "${beds[@]}"; do
        awk -F'\t' '{ if (NF!=6 || $2=="NA" || $3=="NA" || $2+0>=$3+0 || $1=="") exit 1 }' "$b" \
            || { fail "malformed BED: $(basename "$b")"; bad=1; }
    done
    (( bad == 0 )) && pass "every BED row has 6 columns, no NA, start < end"
fi
shopt -u nullglob

# ---------------------------------------------------------------------------
group "Report and manifest"
# ---------------------------------------------------------------------------
assert_min_size "report.html is non-trivial" "$OUT/report.html" 50000
"${RUN[@]}" python - "$OUT/run_manifest.json" <<'PY' && pass "manifest has config, seeds, versions and input checksums" || fail "manifest incomplete"
import json, sys
m = json.load(open(sys.argv[1]))
assert {"config", "seeds", "versions", "inputs", "session_info"} <= set(m), sorted(m)
assert m["config"]["differential"]["fdr"] is not None
assert any(v for v in m["seeds"].values()), "no seeds recorded"
assert all(i["sha256"] for i in m["inputs"]), "input checksum missing"
PY

# ---------------------------------------------------------------------------
group "Golden tables"
# ---------------------------------------------------------------------------
mkdir -p "$GOLDEN"
for arm in $ARMS; do
    for f in "${arm}_supercluster_sizes.tsv" "${arm}_k_diagnostics.tsv"; do
        src="$OUT/clusters/$f"
        gold="$GOLDEN/$f"
        if [[ "${UPDATE_GOLDEN:-0}" == "1" ]]; then
            cp "$src" "$gold"; pass "golden updated: $f"
        elif [[ ! -f "$gold" ]]; then
            cp "$src" "$gold"; warn "golden created (nothing to compare against yet): $f"
        elif "${RUN[@]}" python tests/helpers/compare_tables.py "$gold" "$src" 2>"$SANDBOX/cmp.txt"; then
            # Compared as data, not bytes: floats may differ in the last
            # digits across CPUs and BLAS builds. See compare_tables.py.
            pass "golden matches: $f"
        else
            fail "golden differs: $f. If intended, rerun with UPDATE_GOLDEN=1"
            head -10 "$SANDBOX/cmp.txt"
        fi
    done
done

# ---------------------------------------------------------------------------
group "Idempotency"
# ---------------------------------------------------------------------------
if "${RUN[@]}" snakemake --configfile "$CFG" -j 4 2>&1 | grep -q "Nothing to be done"; then
    pass "second run is a no-op"
else
    fail "second run did work; something is not declaring its outputs correctly"
fi

summary

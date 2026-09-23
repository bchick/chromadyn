#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Tier 1: static checks. Seconds, no pipeline execution, no data required
# beyond the committed demo files.
#
#   bash tests/test_static.sh
#
# Exit status is the result: 0 if nothing failed.
# ---------------------------------------------------------------------------

source "$(dirname "${BASH_SOURCE[0]}")/lib/assert.sh"
cd "$REPO_ROOT" || exit 1

# Run R and snakemake through pixi when we are not already inside the env.
if [[ -n "${CONDA_PREFIX:-}" && -x "${CONDA_PREFIX}/bin/Rscript" ]]; then
    RUN=(); else RUN=(pixi run --frozen)
fi
r()    { "${RUN[@]}" Rscript "$@"; }
snake(){ "${RUN[@]}" snakemake "$@"; }
py()   { "${RUN[@]}" python "$@"; }

CONFIGS=(config/config.yaml config/demo.yaml)
while IFS= read -r f; do CONFIGS+=("$f"); done < <(find config/presets -name '*.yaml' 2>/dev/null | sort)

# ---------------------------------------------------------------------------
group "Workflow parses"
# ---------------------------------------------------------------------------
# Dry runs go into an empty output directory. Against an existing results/,
# Snakemake does not schedule jobs whose outputs are present and so never
# checks their inputs: a missing input can pass here and fail on a clean CI
# runner, which is exactly what happened once.
DRY="$(mktemp -d "${TMPDIR:-/tmp}/timecourse-patterns-dry-XXXXXX")"
trap 'rm -rf "$DRY"' EXIT
printf "output:\n  dir: %s/demo\n" "$DRY" > "$DRY/demo_out.yaml"
py tests/helpers/merge_config.py config/demo.yaml "$DRY/demo_out.yaml" > "$DRY/demo.yaml"
# Gene mode as a real config file. Passing it with --config does not work:
# the value is not parsed as YAML there, so `null` arrives as the string
# "null" and the workflow goes looking for a BED file of that name.
printf "input:\n  features: null\noutput:\n  dir: %s/gene\n" "$DRY" > "$DRY/gene_out.yaml"
py tests/helpers/merge_config.py config/demo.yaml "$DRY/gene_out.yaml" > "$DRY/gene.yaml"

assert_ok "snakemake --lint" snake --lint --configfile "$DRY/demo.yaml"
assert_ok "dry run, demo, into an empty directory" snake -n --configfile "$DRY/demo.yaml"
assert_ok "dry run, gene mode, into an empty directory" snake -n --configfile "$DRY/gene.yaml"

# ---------------------------------------------------------------------------
group "Every script: path exists"
# ---------------------------------------------------------------------------
missing=0
while IFS= read -r s; do
    if [[ -f "workflow/$s" || -f "workflow/rules/$s" ]]; then :; else
        # script: paths in rules/*.smk are relative to the rule file
        resolved="workflow/rules/$s"
        [[ -f "$(cd "$(dirname "$resolved")" 2>/dev/null && pwd)/$(basename "$s")" ]] || missing=$((missing+1))
        [[ -f "$resolved" ]] || fail "script not found: $s"
    fi
done < <(grep -rhoP '(?<=script:)\s*"\K[^"]+' workflow/ | sort -u)
(( missing == 0 )) && pass "all script: paths resolve"

# ---------------------------------------------------------------------------
group "Config and schema"
# ---------------------------------------------------------------------------
for c in "${CONFIGS[@]}"; do
    assert_ok "schema validates $(basename "$c")" py tests/helpers/validate_config.py "$c"
done

# The shipped config must declare every key the schema requires, or the
# documentation and the validation have drifted apart.
assert_ok "config.yaml declares every key the schema requires" \
    py tests/helpers/check_schema_coverage.py

# Paths must resolve identically on the Python side and the R side, since both
# read workflow/paths.yaml and a divergence would be invisible until runtime.
assert_ok "paths.yaml resolves identically in Python and R" \
    "${RUN[@]}" bash tests/helpers/check_paths_agree.sh

# ---------------------------------------------------------------------------
group "R library"
# ---------------------------------------------------------------------------
assert_ok "tests/test_lib.R (guardrail unit tests)" r tests/test_lib.R
for f in workflow/scripts/*.R workflow/scripts/lib/*.R demo/*.R; do
    [[ -e "$f" ]] || continue
    assert_ok "parses: $f" r -e "invisible(parse('$f'))"
done

# ---------------------------------------------------------------------------
group "No hard-coded thresholds outside the config"
# ---------------------------------------------------------------------------
# The single most-violated rule in the project this method came from, and the
# main reason it was not reusable. lib/ is exempt: it holds the fallbacks that
# cfg() defaults to. Comments are stripped before matching so that prose
# mentioning a value is not a violation.
viol=0
for f in workflow/scripts/[0-9]*.R; do
    [[ -e "$f" ]] || continue
    hits=$(sed 's/#.*//' "$f" | grep -nE '(^|[^._[:alnum:]])(0\.01|0\.05|30000|12000|5000)([^._[:alnum:]]|$)' || true)
    if [[ -n "$hits" ]]; then
        while IFS= read -r h; do fail "hard-coded threshold in $f: $h"; viol=$((viol+1)); done <<< "$hits"
    fi
done
(( viol == 0 )) && pass "no hard-coded thresholds in workflow/scripts/*.R"

# Every cfg() key a script reads must exist in config.yaml, or it silently
# takes a default that is documented nowhere.
assert_ok "every cfg() key used by a script exists in config.yaml" \
    py tests/helpers/check_cfg_keys.py

# ---------------------------------------------------------------------------
group "Prose conventions"
# ---------------------------------------------------------------------------
# The pattern is written as a PCRE escape so that this file does not itself
# contain the character it is looking for. docs/provenance/ is exempt: those
# files are verbatim copies of the notebooks this method came from and have to
# stay byte-identical to their source to be worth keeping.
emdash=$(grep -rlP '\x{2014}' --include='*.R' --include='*.smk' --include='*.yaml' \
    --include='*.md' --include='*.sh' --include='*.toml' --include='Snakefile' \
    config workflow tests demo docs README.md 2>/dev/null \
    | grep -vE 'tests/test_static.sh|docs/provenance/' || true)
if [[ -z "$emdash" ]]; then pass "no em dashes"
else while IFS= read -r f; do fail "em dash in $f"; done <<< "$emdash"; fi

ai=$(grep -rlniE 'co-authored-by: *claude|generated with \[claude|🤖 generated' \
    --include='*.R' --include='*.smk' --include='*.yaml' --include='*.md' \
    --include='*.sh' --include='*.cff' . 2>/dev/null \
    | grep -vE 'tests/test_static.sh' || true)
if [[ -z "$ai" ]]; then pass "no AI attribution in tracked files"
else while IFS= read -r f; do fail "AI attribution in $f"; done <<< "$ai"; fi

# ---------------------------------------------------------------------------
group "Optional linters"
# ---------------------------------------------------------------------------
if "${RUN[@]}" shellcheck --version >/dev/null 2>&1; then
    for f in tests/*.sh tests/lib/*.sh tests/helpers/*.sh workflow/envs/*.sh; do
        [[ -e "$f" ]] || continue
        assert_ok "shellcheck $f" "${RUN[@]}" shellcheck -S warning "$f"
    done
else
    skip "shellcheck not installed"
fi

if r -e 'quit(status = !requireNamespace("lintr", quietly = TRUE))' 2>/dev/null; then
    assert_ok "lintr over workflow/scripts" r tests/helpers/run_lintr.R
else
    skip "lintr not installed"
fi

# ---------------------------------------------------------------------------
group "Repository hygiene"
# ---------------------------------------------------------------------------
# Rplots.pdf appears whenever R draws with no device open and is pure noise
# in a diff. It got committed once already.
if git ls-files --error-unmatch Rplots.pdf >/dev/null 2>&1; then
    fail "Rplots.pdf is tracked; it is an artifact, not a file"
else
    pass "no stray Rplots.pdf"
fi

summary

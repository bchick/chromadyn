# ---------------------------------------------------------------------------
# assert.sh: shared assertion helpers for the chromadyn test tiers.
#
# Source this, do not execute it. It sets up:
#   - `set -uo pipefail`, deliberately WITHOUT -e, so a failing assertion is
#     recorded and the suite carries on to find the rest
#   - REPO_ROOT, located from this file rather than from the caller's cwd
#   - TTY-conditional colour
#   - the pass/fail/warn/skip/group counter quartet
#   - summary(), which re-lists every failure and returns the suite's status
#
# Lifted from mcf7_project/tests/test_pipelines.sh, which is where the idiom
# comes from.
#
# End every suite with:  summary; exit $?
# ---------------------------------------------------------------------------
# shellcheck shell=bash
# REPO_ROOT and the colour variables are consumed by the suites that source
# this file, so shellcheck cannot see their use from here.
# shellcheck disable=SC2034

set -uo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

if [[ -t 1 ]]; then
    RED=$'\e[31m'; GREEN=$'\e[32m'; YELLOW=$'\e[33m'; CYAN=$'\e[36m'
    BOLD=$'\e[1m'; DIM=$'\e[2m'; RESET=$'\e[0m'
else
    RED=""; GREEN=""; YELLOW=""; CYAN=""; BOLD=""; DIM=""; RESET=""
fi

PASS=0; FAIL=0; WARN=0; SKIP=0
declare -a FAILED_LINES=()
declare -a WARN_LINES=()

pass()  { PASS=$((PASS+1)); printf "  %sPASS%s %s\n" "$GREEN" "$RESET" "$1"; }
fail()  { FAIL=$((FAIL+1)); printf "  %sFAIL%s %s\n" "$RED" "$RESET" "$1"; FAILED_LINES+=("$1"); }
warn()  { WARN=$((WARN+1)); printf "  %sWARN%s %s\n" "$YELLOW" "$RESET" "$1"; WARN_LINES+=("$1"); }
skip()  { SKIP=$((SKIP+1)); printf "  %sSKIP%s %s\n" "$DIM" "$RESET" "$1"; }
group() { printf "\n%s%s%s\n" "$BOLD" "$1" "$RESET"; }

# assert_ok "label" command...     command must exit 0
assert_ok() {
    local label="$1"; shift
    if "$@" >/dev/null 2>&1; then pass "$label"; else fail "$label"; fi
}

# assert_fails "label" command...  command must exit non-zero
assert_fails() {
    local label="$1"; shift
    if "$@" >/dev/null 2>&1; then fail "$label (expected non-zero exit)"; else pass "$label"; fi
}

# assert_file "label" path         file must exist and be non-empty
assert_file() {
    local label="$1" path="$2"
    if [[ -s "$path" ]]; then pass "$label"
    elif [[ -e "$path" ]]; then fail "$label (exists but empty: $path)"
    else fail "$label (missing: $path)"; fi
}

# assert_min_size "label" path bytes
assert_min_size() {
    local label="$1" path="$2" min="$3"
    if [[ ! -e "$path" ]]; then fail "$label (missing: $path)"; return; fi
    local sz; sz=$(wc -c < "$path")
    if (( sz >= min )); then pass "$label ($sz bytes)"
    else fail "$label ($sz bytes, expected at least $min)"; fi
}

# assert_eq "label" actual expected
assert_eq() {
    local label="$1" actual="$2" expected="$3"
    if [[ "$actual" == "$expected" ]]; then pass "$label"
    else fail "$label (got '$actual', expected '$expected')"; fi
}

summary() {
    printf "\n%s%s%s\n" "$BOLD" "----------------------------------------" "$RESET"
    printf "%sPASS%s %-4d  %sFAIL%s %-4d  %sWARN%s %-4d  %sSKIP%s %d\n" \
        "$GREEN" "$RESET" "$PASS" "$RED" "$RESET" "$FAIL" \
        "$YELLOW" "$RESET" "$WARN" "$DIM" "$RESET" "$SKIP"
    if (( ${#WARN_LINES[@]} )); then
        printf "\n%sWarnings:%s\n" "$YELLOW" "$RESET"
        printf "  - %s\n" "${WARN_LINES[@]}"
    fi
    if (( ${#FAILED_LINES[@]} )); then
        printf "\n%sFailures:%s\n" "$RED" "$RESET"
        printf "  - %s\n" "${FAILED_LINES[@]}"
    fi
    [[ $FAIL -eq 0 ]]
}

#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Finish the Bioconductor install that pixi deliberately leaves undone.
#
# Several bioconda bioconductor-* packages ship their payload in a conda
# post-link script rather than in the package itself. GenomeInfoDbData is the
# one that matters here: GenomeInfoDb depends on it, and so transitively do
# DESeq2, GenomicRanges, rtracklayer and DEGreport. pixi skips post-link
# scripts by default (they are arbitrary code running at install time), so
# without this step `library(DEGreport)` fails with
# "there is no package called 'GenomeInfoDbData'".
#
# Rather than ask every user to flip pixi's run-post-link-scripts setting to
# `insecure`, timecourse-patterns runs the specific scripts it needs, here, where they
# are visible and auditable. The pixi tasks that need R depend on this one.
#
# Idempotent: exits immediately once the packages load.
# ---------------------------------------------------------------------------
set -uo pipefail

if [[ -z "${CONDA_PREFIX:-}" ]]; then
  echo "postinstall: CONDA_PREFIX is not set. Run this through pixi, for example" >&2
  echo "             pixi run postinstall" >&2
  exit 1
fi

# Data packages installed by post-link scripts, as `r_package_name:script_stem`.
NEEDED=(
  "GenomeInfoDbData:.bioconductor-genomeinfodbdata-post-link.sh"
)

missing=()
for entry in "${NEEDED[@]}"; do
  pkg="${entry%%:*}"
  if ! Rscript -e "quit(status = !requireNamespace('${pkg}', quietly = TRUE))" 2>/dev/null; then
    missing+=("$entry")
  fi
done

if [[ ${#missing[@]} -eq 0 ]]; then
  echo "postinstall: Bioconductor data packages already present, nothing to do."
  exit 0
fi

for entry in "${missing[@]}"; do
  pkg="${entry%%:*}"
  script="${entry##*:}"
  path="${CONDA_PREFIX}/bin/${script}"
  if [[ ! -f "$path" ]]; then
    echo "postinstall: expected post-link script not found: $path" >&2
    echo "             The environment may be incomplete; try 'pixi install'." >&2
    exit 1
  fi
  echo "postinstall: installing ${pkg} ..."
  # The script reads $PREFIX, which conda sets at install time and pixi does not.
  PREFIX="$CONDA_PREFIX" bash "$path" >/dev/null 2>&1
  if ! Rscript -e "quit(status = !requireNamespace('${pkg}', quietly = TRUE))" 2>/dev/null; then
    echo "postinstall: ${pkg} still not installed after running ${script}." >&2
    echo "             Re-run without output suppression to see why:" >&2
    echo "               PREFIX=\"\$CONDA_PREFIX\" bash \"$path\"" >&2
    exit 1
  fi
  echo "postinstall: ${pkg} OK"
done

#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Fetch the processed counts for GEO series GSE227634 and check them.
#
#   bash vignettes/data/download_gse227634.sh [outdir]
#
# GSE227634 is the RNA-seq arm of the SuperSeries the bundled demo's ATAC data
# comes from: McDonald, Chick et al., Immunity 56(6):1303-1319.e5 (2023),
# doi:10.1016/j.immuni.2023.05.005.
#
# Downloads are pinned by SHA-256, verified on 2026-09-22. A checksum mismatch
# means GEO reprocessed the deposit; do not silently accept it.
# ---------------------------------------------------------------------------
set -uo pipefail

OUT="${1:-vignettes/data/gse227634}"
BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE227nnn/GSE227634/suppl"

FILES=(
  "GSE227634_featureCounts.d3.mm39.vM27.txt.gz a5120f5bf1f32faf3e2d47539778aac369512879724c1b1150f230cdd51e390a"
  "GSE227634_featureCounts.d5.mm39.vM27.WT.txt.gz 7ae388662dca04565519b049a5917d719e5ae95e47afb7daa8f43ff95a2e4c15"
  "GSE227634_featureCounts.d8.mm39.vM27.final.txt.gz 62f15a1c20da1c8869b416b8590c41d4c196abb33be5bc499b310002119a01ac"
)

mkdir -p "$OUT"
status=0
for entry in "${FILES[@]}"; do
    name="${entry%% *}"; want="${entry##* }"
    dest="$OUT/$name"
    if [[ ! -f "$dest" ]]; then
        echo "downloading $name"
        curl -fsSL -o "$dest" "$BASE/$name" || { echo "  download failed" >&2; status=1; continue; }
    fi
    got=$(sha256sum "$dest" | cut -d' ' -f1)
    if [[ "$got" == "$want" ]]; then
        echo "  ok  $name"
    else
        echo "  CHECKSUM MISMATCH $name" >&2
        echo "    expected $want" >&2
        echo "    got      $got" >&2
        status=1
    fi
done
exit $status

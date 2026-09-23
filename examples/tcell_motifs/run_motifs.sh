#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# run_motifs.sh: known-motif enrichment per chromadyn trajectory class.
#
#   bash examples/tcell_motifs/run_motifs.sh [results_dir] [genome.fa] [motifs.meme]
#
# Each class is tested with MEME-suite SEA against a background of STATIC
# peaks: accessible regions that the same run judged not to change over time.
# That background matters. Against shuffled or genomic sequence, every class
# would be enriched for the motifs of anything that makes chromatin
# accessible; against static peaks, enrichment means "specific to this
# temporal program".
#
# All regions are trimmed to 200 bp around the peak centre, so that peak
# width, which differs between classes, cannot masquerade as enrichment.
#
# Needs: bedtools, MEME suite (sea), a genome FASTA with a .fai index, and a
# MEME-format motif file (JASPAR2024 CORE vertebrates non-redundant was used;
# see README.md for the download and its checksum).
# ---------------------------------------------------------------------------
set -uo pipefail

RES="${1:-examples/tcell_motifs/work/results}"
FA="${2:-/data/bchick/wproj/tcell_project/data/reference/GRCm39.primary_assembly.genome.fa}"
MEME="${3:-examples/tcell_motifs/work/motifs/jaspar2024_vert_nr.meme}"
ARM="WT"
HALF=100
N_BG=20000
SEED=42
OUT="$(dirname "$RES")/motifs"
mkdir -p "$OUT"

for f in "$FA" "$FA.fai" "$MEME" "$RES/clusters/${ARM}_clusters.tsv" \
         "$RES/differential/${ARM}_results.tsv.gz"; do
    [[ -f "$f" ]] || { echo "missing: $f" >&2; exit 1; }
done

centre() {   # BED6 on stdin -> 200 bp window around the centre
    awk -v h=$HALF 'BEGIN{OFS="\t"} {c=int(($2+$3)/2); s=c-h; if (s<0) s=0; print $1,s,c+h,$4}'
}

# --- background: static peaks, subsampled with a fixed seed -----------------
zcat "$RES/differential/${ARM}_results.tsv.gz" | awk -F'\t' 'NR>1 && $NF=="FALSE" {print $1}' \
    > "$OUT/static_ids.txt"
awk 'NR==FNR {k[$1]; next} ($4 in k)' "$OUT/static_ids.txt" \
    examples/tcell_motifs/work/inputs/features.bed \
    | shuf -n $N_BG --random-source=<(yes $SEED) | centre > "$OUT/background.bed"
bedtools getfasta -fi "$FA" -bed "$OUT/background.bed" -fo "$OUT/background.fa"
echo "background: $(grep -c '^>' "$OUT/background.fa") static peaks"

# --- one SEA run per class -------------------------------------------------
tail -n +2 "$RES/clusters/${ARM}_clusters.tsv" | awk -F'\t' '$4!="Unassigned"' \
    | cut -f4 | sort -u | while IFS= read -r cls; do
    slug=$(echo "$cls" | tr 'A-Z ' 'a-z_')
    tail -n +2 "$RES/clusters/${ARM}_clusters.tsv" \
        | awk -F'\t' -v c="$cls" 'BEGIN{OFS="\t"} $4==c {print $7,$8,$9,$1,0,"."}' \
        | centre > "$OUT/${slug}.bed"
    bedtools getfasta -fi "$FA" -bed "$OUT/${slug}.bed" -fo "$OUT/${slug}.fa"
    sea --p "$OUT/${slug}.fa" --n "$OUT/background.fa" --m "$MEME" \
        --o "$OUT/sea_${slug}" --seed $SEED >/dev/null 2>&1 \
        || { echo "sea failed for $cls" >&2; exit 1; }
    printf "  %-22s %6d peaks  %4d enriched motifs\n" "$cls" \
        "$(grep -c '^>' "$OUT/${slug}.fa")" \
        "$(grep -vc '^#\|^RANK\|^$' "$OUT/sea_${slug}/sea.tsv")"
done

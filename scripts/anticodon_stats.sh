#!/usr/bin/env bash
# Per-sample anticodon-count table from tRNAscan-SE outputs (cytosolic + mito).
# Usage: anticodon_stats.sh <sample> <cyto_trnascan.out> <mito_trnascan.out> [outdir]
set -euo pipefail
SAMPLE=$1; CYTO=$2; MITO=$3; OUTDIR=${4:-.}
TSV="$OUTDIR/${SAMPLE}.anticodon_counts.tsv"
{
  printf "compartment\tAA\tanticodon\tcount\n"
  tail -n +4 "$CYTO" | awk -F'\t' '{a=$5;c=$6;gsub(/ /,"",a);gsub(/ /,"",c); if(c!="")print a"\t"c}' \
    | sort | uniq -c | awk '{print "cyto\t"$2"\t"$3"\t"$1}' | sort -k2,2 -k4,4nr
  tail -n +4 "$MITO" | awk -F'\t' '{a=$5;c=$6;gsub(/ /,"",a);gsub(/ /,"",c); if(c!="")print a"\t"c}' \
    | sort -u | awk '{print "mito\t"$1"\t"$2"\t1"}' | sort -k2,2
} > "$TSV"
echo "wrote $TSV"
echo "## $SAMPLE  cytosolic anticodon counts"
awk -F'\t' '$1=="cyto"{print}' "$TSV" | awk -F'\t' '{printf "%s-%s:%s  ",$2,$3,$4}'; echo
echo "cyto: $(awk -F"\t" '$1=="cyto"{n+=$4}END{print n}' "$TSV") tRNAs across $(awk -F"\t" '$1=="cyto"' "$TSV"|wc -l) anticodons / $(awk -F"\t" '$1=="cyto"{print $2}' "$TSV"|sort -u|wc -l) amino acids"
echo "mito: $(awk -F"\t" '$1=="mito"' "$TSV"|wc -l) distinct isotypes"

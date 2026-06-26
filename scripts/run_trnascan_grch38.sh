#!/usr/bin/env bash
# Stage 1 — build the GRCh38 reference tRNA set by scanning the genome FASTA
# with tRNAscan-SE 2.0 (eukaryotic mode, Infernal + tRNA covariance models).
set -euo pipefail

ROOT=/mnt/data/igor-shared/trnascan
ENV="$ROOT/env"
REF="${GENOME:?set \$GENOME to your GRCh38 genome FASTA}"   # optional: GtRNAdb already publishes this scan
OUT="${OUT:-$ROOT/grch38_ref}"
THREADS="${THREADS:-32}"

export PATH="$ENV/bin:$PATH"
mkdir -p "$OUT"
cd "$OUT"

echo "[$(date)] START tRNAscan-SE on $REF  (threads=$THREADS)"
tRNAscan-SE -E -Q -d \
  --thread "$THREADS" \
  -o "$OUT/grch38.trnascan.out" \
  -f "$OUT/grch38.trnascan.ss" \
  -b "$OUT/grch38.trna.bed" \
  -j "$OUT/grch38.trna.gff3" \
  -a "$OUT/grch38.trna.fa" \
  -m "$OUT/grch38.trnascan.stats" \
  -l "$OUT/grch38.trnascan.log" \
  "$REF"
echo "[$(date)] DONE"
echo -n "predicted tRNA sequences: "; grep -c '^>' "$OUT/grch38.trna.fa" || true

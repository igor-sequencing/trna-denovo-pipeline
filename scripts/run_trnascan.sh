#!/usr/bin/env bash
# General tRNAscan-SE runner. MODE selects search mode, e.g. "-E" or "-M mammal".
set -euo pipefail
ROOT=/mnt/data/igor-shared/trnascan
ENV="$ROOT/env"
REF="${GENOME:?set \$GENOME to the FASTA to scan}"
MODE="${MODE:--E}"
OUT="${OUT:?set OUT}"
PREFIX="${PREFIX:-grch38}"
THREADS="${THREADS:-32}"
export PATH="$ENV/bin:$PATH"
mkdir -p "$OUT"; cd "$OUT"
echo "[$(date)] START tRNAscan-SE $MODE on $REF (threads=$THREADS) -> $OUT"
# shellcheck disable=SC2086
tRNAscan-SE $MODE -Q -d --thread "$THREADS" \
  -o "$OUT/$PREFIX.trnascan.out" \
  -f "$OUT/$PREFIX.trnascan.ss" \
  -b "$OUT/$PREFIX.trna.bed" \
  -j "$OUT/$PREFIX.trna.gff3" \
  -a "$OUT/$PREFIX.trna.fa" \
  -m "$OUT/$PREFIX.trnascan.stats" \
  -l "$OUT/$PREFIX.trnascan.log" \
  "$REF"
echo "[$(date)] DONE $MODE"
echo -n "tRNA sequences: "; grep -c '^>' "$OUT/$PREFIX.trna.fa" 2>/dev/null || echo 0

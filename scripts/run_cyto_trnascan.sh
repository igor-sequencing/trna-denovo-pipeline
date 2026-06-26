#!/usr/bin/env bash
set -euo pipefail
trap 'echo "CYTO EXIT $?"' EXIT
export PATH=/mnt/data/igor-shared/trnascan/env/bin:$PATH
SRC=/mnt/data/igor-shared/trnascan/saute_SQ326CT3_k31/saute_all_variants.fa
OUT=/mnt/data/igor-shared/trnascan/cyto_trnascan_SQ326CT3
mkdir -p "$OUT"; cd "$OUT"
seqkit grep -v -n -r -p 'MT\.trna' "$SRC" > cyto_variants.fa
echo "cytosolic variants in: $(grep -c '^>' cyto_variants.fa)"
rm -f cyto_trnascan.*
tRNAscan-SE -E -Q --thread 8 -o cyto_trnascan.out -f cyto_trnascan.ss -a cyto_trnascan.fa -b cyto_trnascan.bed cyto_variants.fa
echo "[$(date)] CYTO tRNAscan done"
echo "tRNAs called: $(tail -n +4 cyto_trnascan.out | grep -c .)"

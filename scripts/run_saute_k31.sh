#!/usr/bin/env bash
set -euo pipefail
trap 'echo "SAUTE2 EXIT $?"' EXIT
export PATH=/mnt/data/igor-shared/trnascan/env/bin:$PATH
OUT=/mnt/data/igor-shared/trnascan/saute_SQ326CT3_k31
RD=/mnt/data/igor-shared/trnascan/retrieve_SQ326CT3
REF=/mnt/data/igor-shared/trnascan/ref_build/human_trna_ref.fa
mkdir -p "$OUT"; cd "$OUT"
CORES=${CORES:-12}; KMER=${KMER:-31}; SKMER=${SKMER:-21}; MAXV=${MAXV:-1000000}
echo "[$(date)] SAUTE START kmer=$KMER secondary_kmer=$SKMER max_variants=$MAXV cores=$CORES"
saute --cores "$CORES" --kmer "$KMER" --secondary_kmer "$SKMER" --max_variants "$MAXV" \
  --targets "$REF" \
  --reads "$RD/SQ326CT3_trna_R1.fq.gz,$RD/SQ326CT3_trna_R2.fq.gz" \
  --all_variants "$OUT/saute_all_variants.fa" \
  --selected_variants "$OUT/saute_selected_variants.fa" \
  --gfa "$OUT/saute.gfa"
echo "[$(date)] SAUTE FINISHED"
echo "all_variants: $(grep -c '^>' saute_all_variants.fa 2>/dev/null||echo 0)  selected: $(grep -c '^>' saute_selected_variants.fa 2>/dev/null||echo 0)"
echo "mito-derived (MT.trna): $(grep -c 'MT.trna' saute_all_variants.fa 2>/dev/null||echo 0)"

#!/usr/bin/env bash
# Retrieve tRNA-derived read pairs from WGS FASTQ via BBDuk k-mer match
# against the combined human (cytosolic+mito) tRNA reference.
set -euo pipefail
export PATH=/mnt/data/igor-shared/trnascan/env/bin:$PATH
RB=/mnt/data/igor-shared/trnascan/ref_build
OUT=/mnt/data/igor-shared/trnascan/retrieve_SQ326CT3
mkdir -p "$OUT"; cd "$OUT"
R1=${R1:-/mnt/data/igor-shared/LT/SQ326CT3_1.fq.gz}
R2=${R2:-/mnt/data/igor-shared/LT/SQ326CT3_2.fq.gz}
REF=${REF:-$RB/human_trna_ref.fa}
K=${K:-25}
THREADS=${THREADS:-24}

echo "[$(date)] BBDuk tRNA-read retrieval START  k=$K threads=$THREADS"
echo "  R1=$R1"; echo "  R2=$R2"; echo "  ref=$REF"
bbduk.sh -Xmx16g threads="$THREADS" \
  in1="$R1" in2="$R2" \
  ref="$REF" k="$K" rcomp=t \
  outm1="$OUT/SQ326CT3_trna_R1.fq.gz" outm2="$OUT/SQ326CT3_trna_R2.fq.gz" \
  stats="$OUT/bbduk_stats.txt"
echo "[$(date)] BBDuk DONE"
echo "=== retrieved read counts ==="
seqkit stats "$OUT/SQ326CT3_trna_R1.fq.gz" "$OUT/SQ326CT3_trna_R2.fq.gz"

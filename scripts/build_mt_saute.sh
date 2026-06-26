#!/usr/bin/env bash
set -euo pipefail
trap 'echo "MTSAUTE EXIT $?"' EXIT
export PATH=/mnt/data/igor-shared/trnascan/env/bin:$PATH
OUT=/mnt/data/igor-shared/trnascan/mt_saute_SQ326CT3
RB=/mnt/data/igor-shared/trnascan/ref_build
BAM=/mnt/data/igor-shared/IT/SQ326CT3.bam
MTREF=$RB/mt.fa
CORES=${CORES:-8}
mkdir -p "$OUT"; cd "$OUT"

# 1) detect MT contig (pipefail off: awk 'exit' SIGPIPEs samtools, which is harmless here)
set +o pipefail
MTNAME=$(samtools view -H "$BAM" | awk '/^@SQ/{for(i=1;i<=NF;i++) if($i ~ /^SN:(chrMT|chrM|MT|M)$/){print substr($i,4); exit}}')
set -o pipefail
[ -z "${MTNAME:-}" ] && { echo "ERROR: no MT contig in BAM header"; exit 2; }
echo "[$(date)] MT contig in BAM = '$MTNAME'"

# 2) extract ALL MT reads -> mate-sorted paired fastq
echo "[$(date)] extracting MT reads"
samtools view -b -@ "$CORES" "$BAM" "$MTNAME" > mt_reads.bam
samtools collate -@ "$CORES" -O -u mt_reads.bam | \
  samtools fastq -@ "$CORES" -1 mt_R1.fq.gz -2 mt_R2.fq.gz -s mt_singletons.fq.gz -0 /dev/null -n
echo "MT read pairs (R1 lines/4): $(zcat mt_R1.fq.gz | wc -l | awk '{print $1/4}')"

# 3) assemble whole MT with SAUTE
echo "[$(date)] SAUTE assembling whole MT (target=$MTREF)"
saute --cores "$CORES" \
  --targets "$MTREF" \
  --reads mt_R1.fq.gz,mt_R2.fq.gz \
  --all_variants mt_all_variants.fa \
  --selected_variants mt_assembly.fa \
  --gfa mt.gfa
echo "MT assembly:"; seqkit stats mt_assembly.fa

# 4) tRNAscan-SE -M mammal on assembled MT
echo "[$(date)] tRNAscan-SE -M mammal on assembled MT"
rm -f mt_asm_trnascan.*
tRNAscan-SE -M mammal -Q -a mt_asm_trnascan.fa -o mt_asm_trnascan.out -f mt_asm_trnascan.ss mt_assembly.fa
echo "[$(date)] DONE"
echo "mito tRNAs in assembly: $(grep -c '^>' mt_asm_trnascan.fa 2>/dev/null || echo 0)"

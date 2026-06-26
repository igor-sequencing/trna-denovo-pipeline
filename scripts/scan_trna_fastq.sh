#!/usr/bin/env bash
# Lean FASTQ-native tRNA anticodon-switch pipeline (no whole-FASTQ alignment, no SAUTE):
#   FASTQ --BBDuk(paired; either-mate match keeps the pair)--> tRNA reads
#         --bwa-mem2 of ONLY those reads to the full genome--> small indexed BAM
#         --pileup anticodon positions (scan_noncanonical_anticodons.sh)--> switch calls
# Genome (GRCh38 + bwa-mem2 index) supplied via $GENOME. tRNA ref/.ss from $RB.
set -euo pipefail
trap 'echo "FASTQ_PIPE EXIT $?"' EXIT
export PATH=/mnt/data/igor-shared/trnascan/env/bin:$PATH
R1=${1:?usage: scan_trna_fastq.sh <R1.fq.gz> <R2.fq.gz>}; R2=${2:?need R2}
GENOME=${GENOME:?set \$GENOME to a GRCh38 FASTA with a bwa-mem2 index}
RB=/mnt/data/igor-shared/trnascan/ref_build
TRNAREF=$RB/human_trna_ref.fa
CORES=${CORES:-16}; K=${K:-23}
SAMPLE=${SAMPLE:-$(basename "$R1" | sed "s/[._]R\{0,1\}1.*//")}_fq
OUT=/mnt/data/igor-shared/trnascan/fastq_${SAMPLE}
mkdir -p "$OUT"; cd "$OUT"

echo "[$(date)] STEP1 BBDuk k-mer filter (paired -> mates kept)"
bbduk.sh -Xmx16g threads="$CORES" in1="$R1" in2="$R2" ref="$TRNAREF" k="$K" rcomp=t \
  outm1=trna_R1.fq.gz outm2=trna_R2.fq.gz stats=bbduk_stats.txt 2>bbduk.log
echo "  filtered pairs: $(zcat trna_R1.fq.gz | wc -l | awk '{print $1/4}')"

echo "[$(date)] STEP2 align ONLY the filtered pairs to the full genome"
bwa-mem2 mem -t "$CORES" "$GENOME" trna_R1.fq.gz trna_R2.fq.gz 2>bwa.log | \
  samtools sort -@ "$CORES" -o "$SAMPLE.bam" -
samtools index "$SAMPLE.bam"
echo "  primary reads MAPQ>=20: $(samtools view -c -F0x900 -q20 "$SAMPLE.bam")"

echo "[$(date)] STEP3 anticodon-switch pileup"
/mnt/data/igor-shared/trnascan/scan_noncanonical_anticodons.sh "$OUT/$SAMPLE.bam" 3 10
echo "[$(date)] FASTQ_PIPE DONE"

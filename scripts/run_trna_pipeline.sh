#!/usr/bin/env bash
# End-to-end tRNA de-novo pipeline from a BAM: cytosolic (locus extract -> SAUTE k31 -> tRNAscan -E)
# + mito (whole-MT extract -> SAUTE -> tRNAscan -M mammal) + anticodon lists + non-canonical screen.
set -euo pipefail
trap 'echo "PIPELINE EXIT $?"' EXIT
export PATH=/mnt/data/igor-shared/trnascan/env/bin:$PATH
BAM=${1:?usage: run_trna_pipeline.sh <bam>}
RB=/mnt/data/igor-shared/trnascan/ref_build   # built genome-free by build_ref.sh
CYTOREF=$RB/hg38-tRNAs.fa
MTREF=$RB/mt.fa
CORES=${CORES:-12}
SAMPLE=$(basename "$BAM" .bam)
OUT=/mnt/data/igor-shared/trnascan/pipeline_${SAMPLE}
mkdir -p "$OUT"; cd "$OUT"
echo "[$(date)] PIPELINE start sample=$SAMPLE"

# 0) contig naming + MT name
HDR=$(samtools view -H "$BAM")
printf '%s\n' "$HDR" | awk '/^@SQ/{for(i=1;i<=NF;i++) if($i ~ /^SN:/) print substr($i,4)}' | sort -u > bam_contigs.txt
MTNAME=$(awk '$1=="MT"||$1=="chrMT"||$1=="chrM"||$1=="M"{print;exit}' bam_contigs.txt)
[ -z "${MTNAME:-}" ] && { echo "ERROR no MT contig"; exit 2; }
if grep -qE '^chr[0-9XY]' bam_contigs.txt; then CHRPFX=yes; else CHRPFX=no; fi
echo "MT=$MTNAME chr-prefixed=$CHRPFX"

# 1) cytosolic tRNA loci BED reconciled to BAM naming, kept only on present contigs
if [ "$CHRPFX" = no ]; then sed 's/^chr//' "$RB/hg38-tRNAs.bed" > loci.rec.bed; else cp "$RB/hg38-tRNAs.bed" loci.rec.bed; fi
awk 'NR==FNR{c[$1]=1;next} ($1 in c)' bam_contigs.txt loci.rec.bed > trna_loci.bed
echo "cytosolic loci usable: $(wc -l < trna_loci.bed)/$(wc -l < loci.rec.bed)"

# 2) cytosolic reads (all alignments overlapping loci, primary only) -> single-end fastq
echo "[$(date)] extract cytosolic reads"
samtools view -@ "$CORES" -F 0x900 -L trna_loci.bed "$BAM" | awk '{print "@r"NR"\n"$10"\n+\n"$11}' | gzip > cyto_all.fq.gz
echo "cyto reads: $(zcat cyto_all.fq.gz | wc -l | awk '{print $1/4}')"

# 3) mito reads (all on MT) -> paired fastq
echo "[$(date)] extract mito reads"
samtools view -b -@ "$CORES" "$BAM" "$MTNAME" > mt_reads.bam
samtools collate -@ "$CORES" -O -u mt_reads.bam | samtools fastq -@ "$CORES" -1 mt_R1.fq.gz -2 mt_R2.fq.gz -s /dev/null -0 /dev/null -n 2>/dev/null
echo "mito pairs: $(zcat mt_R1.fq.gz | wc -l | awk '{print $1/4}')"

# 4) SAUTE cytosolic (kmer 31) + mito (auto)
echo "[$(date)] SAUTE cytosolic"
saute --cores "$CORES" --kmer 31 --secondary_kmer 21 --max_variants 1000000 \
  --targets "$CYTOREF" --reads cyto_all.fq.gz \
  --all_variants cyto_variants.fa --selected_variants cyto_selected.fa --gfa cyto.gfa || echo "WARN cyto-saute rc=$?"
echo "[$(date)] SAUTE mito whole-MT"
saute --cores "$CORES" --targets "$MTREF" --reads mt_R1.fq.gz,mt_R2.fq.gz \
  --all_variants mt_all_variants.fa --selected_variants mt_assembly.fa --gfa mt.gfa || echo "WARN mt-saute rc=$?"

# 5) tRNAscan
echo "[$(date)] tRNAscan -E cytosolic"
tRNAscan-SE -E -Q --thread "$CORES" -o cyto_trnascan.out -f cyto_trnascan.ss -a cyto_trnascan.fa cyto_variants.fa || echo "WARN cyto-scan rc=$?"
echo "[$(date)] tRNAscan -M mammal mito"
tRNAscan-SE -M mammal -Q -o mt_trnascan.out -f mt_trnascan.ss -a mt_trnascan.fa mt_assembly.fa || echo "WARN mt-scan rc=$?"

# 6) summary
echo "===== RESULTS $SAMPLE ====="
echo "cyto variants: $(grep -c '^>' cyto_variants.fa 2>/dev/null||echo 0)"
echo "cyto tRNAs called: $(tail -n +4 cyto_trnascan.out 2>/dev/null|grep -c .)"
echo "cyto distinct anticodons: $(tail -n +4 cyto_trnascan.out 2>/dev/null|awk -F'\t' '{gsub(/ /,"",$6);print $6}'|sort -u|grep -c .)"
echo "mito assembly contigs: $(grep -c '^>' mt_assembly.fa 2>/dev/null||echo 0)"
echo "mito tRNA calls: $(grep -c '^>' mt_trnascan.fa 2>/dev/null||echo 0)"
echo "mito distinct isotypes: $(tail -n +4 mt_trnascan.out 2>/dev/null|awk -F'\t' '{gsub(/ /,"",$6);print $6}'|sort -u|grep -c .)"
grep '^>' "$CYTOREF" | sed 's/.*tRNA-//' | awk -F- '{print $2}' | grep -vi NNN | sort -u > ref_ac.txt
tail -n +4 cyto_trnascan.out 2>/dev/null | awk -F'\t' '{gsub(/ /,"",$6);print $6}' | grep -v '^$' | sort -u > asm_ac.txt
echo "NON-CANONICAL cyto anticodons (assembly minus canonical): $(comm -23 asm_ac.txt ref_ac.txt | tr '\n' ' ')"
echo "[$(date)] PIPELINE done"

#!/usr/bin/env bash
# Build the combined human tRNA reference (GENOME-FREE):
#   cytosolic = GtRNAdb hg38 (public download; tRNAscan-SE on GRCh38)
#   mito      = 22 tRNAs from the public rCRS / NC_012920.1 via tRNAscan-SE -M mammal
# Outputs in $RB: human_trna_ref.fa, hg38-tRNAs.fa, hg38-tRNAs-detailed.ss (scanner input),
#                 rCRS.NC_012920.1.fa, mito-tRNAs.fa
set -euo pipefail
export PATH=/mnt/data/igor-shared/trnascan/env/bin:$PATH
RB=${RB:-/mnt/data/igor-shared/trnascan/ref_build}
mkdir -p "$RB"; cd "$RB"

echo "=== cytosolic: GtRNAdb hg38 (public) ==="
curl -sL -o hg38-tRNAs.tar.gz https://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/hg38-tRNAs.tar.gz
tar xzf hg38-tRNAs.tar.gz
echo "cytosolic tRNAs: $(grep -c '^>' hg38-tRNAs.fa)"

echo "=== mitochondrial: public rCRS (NC_012920.1) — no genome reference needed ==="
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_012920.1&rettype=fasta&retmode=text" -o rCRS.NC_012920.1.fa
rm -f mito-tRNAs.fa mito.out mito.ss
tRNAscan-SE -M mammal -Q -a mito-tRNAs.fa -o mito.out -f mito.ss rCRS.NC_012920.1.fa
echo "mito tRNAs: $(grep -c '^>' mito-tRNAs.fa)"

echo "=== combine ==="
cat hg38-tRNAs.fa mito-tRNAs.fa > human_trna_ref.fa
echo "combined: $(grep -c '^>' human_trna_ref.fa)"
seqkit stats human_trna_ref.fa

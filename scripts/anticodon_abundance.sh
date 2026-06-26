#!/usr/bin/env bash
# Per-anticodon read-count abundance (Poisson unit) for an overdispersion test.
# reads_a ~= (sum of bedcov base-depth over the anticodon's loci) / read_length
# Usage: anticodon_abundance.sh <bam> <label>
set -uo pipefail
export PATH=/mnt/data/igor-shared/trnascan/env/bin:$PATH
BAM=$1; LABEL=$2
RB=/mnt/data/igor-shared/trnascan/ref_build
LBED=$RB/trna_loci_labeled.bed
RL=$(samtools view -F0x900 "$BAM" 2>/dev/null | head -1 | awk "{print length(\$10)}")
samtools bedcov "$LBED" "$BAM" | awk -v S="$LABEL" -v rl="$RL" '
{ ds[$4]+=$5 } END{ for(a in ds) printf "%s\t%s\t%.1f\n", S, a, ds[a]/rl }'

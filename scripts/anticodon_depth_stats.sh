#!/usr/bin/env bash
# Depth-based normalized coverage per anticodon (dosage proxy, NOT assembly multiplicity).
# For each anticodon: mean read depth across its tRNA gene loci / sample genome-wide mean depth.
# Usage: anticodon_depth_stats.sh <bam> <label>
set -uo pipefail
export PATH=/mnt/data/igor-shared/trnascan/env/bin:$PATH
BAM=$1; LABEL=$2
RB=/mnt/data/igor-shared/trnascan/ref_build
LBED=$RB/trna_loci_labeled.bed
[ -f "$LBED" ] || awk -F"\t" '{n=$4; sub(/^tRNA-/,"",n); split(n,a,"-"); aa=a[1]; ac=a[2];
  if(ac=="NNN"||ac=="") next; c=$1; sub(/^chr/,"",c); print c"\t"$2"\t"$3"\t"aa"-"ac}' "$RB/hg38-tRNAs.bed" > "$LBED"
RL=$(samtools view -F0x900 "$BAM" 2>/dev/null | head -1 | awk "{print length(\$10)}")
read MAP LEN < <(samtools idxstats "$BAM" | awk "{m+=\$3; l+=\$2} END{print m, l}")
G=$(awk -v m="$MAP" -v rl="$RL" -v l="$LEN" "BEGIN{printf \"%.4f\", m*rl/l}")
samtools bedcov "$LBED" "$BAM" | awk -v G="$G" -v S="$LABEL" "
{ ac=\$4; ds[ac]+=\$5; ln[ac]+=(\$3-\$2) }
END{ for(a in ds){ md=ds[a]/ln[a]; printf \"%s\t%s\t%.2f\t%.3f\n\", S, a, md, md/G } }"

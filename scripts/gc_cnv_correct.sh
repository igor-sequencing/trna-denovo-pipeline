#!/usr/bin/env bash
# GC-corrected per-anticodon copy number for one sample.
# baseline(GC) = mean control-window depth at that GC (genome-wide); corrected_locus = depth/baseline(GC).
set -uo pipefail
export PATH=/mnt/data/igor-shared/trnascan/env/bin:$PATH
BAM=$1; L=$2
WD=/mnt/data/igor-shared/trnascan/gc_cnv; cd "$WD"
RB=/mnt/data/igor-shared/trnascan/ref_build
samtools bedcov control.bed "$BAM" > ctrl_$L.cov
samtools bedcov "$RB/trna_loci_labeled.bed" "$BAM" > trna_$L.cov
awk -v S="$L" '
FILENAME ~ /control_gc/ { gc[$1]=$2; next }
FILENAME ~ /ctrl_.*\.cov/ { k=$1":"$2; d=$4/200; g=gc[k]; if(g!="" && d>0 && d<300){ csum[g]+=d; ccnt[g]++ } next }
FILENAME ~ /trna_gc/   { tg[$1]=$3; tl[$1]=$2; next }
FILENAME ~ /trna_.*\.cov/ { k=$1":"$2; d=$5/($3-$2); g=tg[k]; lab=tl[k];
   if(g!="" && (g in csum) && ccnt[g]>0){ base=csum[g]/ccnt[g]; if(base>0){ tot[lab]+=d/base; n[lab]++ } } next }
END{ for(a in tot) printf "%s\t%s\t%.2f\t%d\n", S, a, tot[a], n[a] }
' control_gc.tsv ctrl_$L.cov trna_gc.tsv trna_$L.cov

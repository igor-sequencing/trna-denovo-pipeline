#!/usr/bin/env bash
# Genome-free non-canonical cytosolic anticodon screen from a BAM.
# Reports EVERY anticodon-position alt allele supported by >= MIN_ALT quality reads
# (MAPQ>=MINMAPQ, baseQ>=20). No VAF floor: real-vs-error is judged by read quality &
# strand balance (reported), not by allele-fraction magnitude. Genome-free (uses public .ss).
# Usage: scan_noncanonical_anticodons.sh <bam> [min_alt_reads=3] [min_locus_depth=10]
set -euo pipefail
export PATH=/mnt/data/igor-shared/trnascan/env/bin:$PATH
BAM=${1:?usage: <bam> [min_alt_reads] [min_locus_depth]}
MINALT=${2:-3}; MINDP=${3:-10}; CORES=${CORES:-8}; MINMAPQ=${MINMAPQ:-20}
RB=/mnt/data/igor-shared/trnascan/ref_build
SS=$RB/hg38-tRNAs-detailed.ss
SAMPLE=$(basename "$BAM" .bam)
OUT=/mnt/data/igor-shared/trnascan/noncanon_scan_${SAMPLE}
mkdir -p "$OUT"; cd "$OUT"

samtools view -H "$BAM" | awk '/^@SQ/{for(i=1;i<=NF;i++) if($i~/^SN:/) print substr($i,4)}' | sort -u > bam_contigs.txt
grep -qE '^chr[0-9XY]' bam_contigs.txt && BAMCHR=yes || BAMCHR=no
awk -v bamchr="$BAMCHR" '
$0 ~ /trna[0-9]+ \(/ && $2 ~ /^\(/ { gene=$1; next }
$0 ~ /Anticodon:/ {
  ac="";coords=""; for(i=1;i<=NF;i++){ if($i=="Anticodon:")ac=$(i+1); if($i~/^\(/)coords=$i }
  gsub(/[()]/,"",coords); split(coords,c,"-"); c1=c[1]+0;c2=c[2]+0;
  if(ac ~ /[^ACGT]/ || c1<1 || c2<1) next;
  lo=(c1<c2)?c1:c2; hi=(c1<c2)?c2:c1; strand=(c1<c2)?"+":"-";
  split(gene,g,"."); ec=g[1]; if(bamchr=="no") sub(/^chr/,"",ec);
  print gene"\t"ec"\t"strand"\t"ac"\t"lo"\t"(lo+1)"\t"hi }' "$SS" > ac_table.tsv
awk 'NR==FNR{c[$1]=1;next} ($2 in c)' bam_contigs.txt ac_table.tsv > ac_table.bam.tsv
awk '{print $2"\t"($5-1)"\t"$5"\n"$2"\t"($6-1)"\t"$6"\n"$2"\t"($7-1)"\t"$7}' ac_table.bam.tsv | sort -k1,1 -k2,2n -u > ac_positions.bed

samtools view -b -@ "$CORES" -q "$MINMAPQ" -L ac_positions.bed "$BAM" > sub.bam
samtools index sub.bam
# per-position strand-split base counts (Q>=20 reads): Af Ar Cf Cr Gf Gr Tf Tr
samtools mpileup -B -q "$MINMAPQ" -Q20 -l ac_positions.bed sub.bam 2>/dev/null | awk '
{ b=$5;n=length(b);split("Af Ar Cf Cr Gf Gr Tf Tr",hdr," ");for(x=1;x<=8;x++)v[x]=0;i=1;
  while(i<=n){ c=substr(b,i,1);
    if(c=="^"){i+=2;continue} else if(c=="$"){i++;continue} else if(c=="*"||c=="#"){i++;continue}
    else if(c=="+"||c=="-"){ j=i+1;num=""; while(substr(b,j,1)~/[0-9]/){num=num substr(b,j,1);j++} i=j+(num+0);continue }
    else { if(c=="A")v[1]++;else if(c=="a")v[2]++;else if(c=="C")v[3]++;else if(c=="c")v[4]++;
           else if(c=="G")v[5]++;else if(c=="g")v[6]++;else if(c=="T")v[7]++;else if(c=="t")v[8]++; i++ } }
  print $1":"$2"\t"v[1]"\t"v[2]"\t"v[3]"\t"v[4]"\t"v[5]"\t"v[6]"\t"v[7]"\t"v[8] }' > pos_counts.tsv

awk -v MINALT="$MINALT" -v MINDP="$MINDP" '
function comp(x){return x=="A"?"T":x=="T"?"A":x=="G"?"C":x=="C"?"G":"N"}
function tot(p){return Af[p]+Ar[p]+Cf[p]+Cr[p]+Gf[p]+Gr[p]+Tf[p]+Tr[p]}
function fwd(p,b){return b=="A"?Af[p]:b=="C"?Cf[p]:b=="G"?Gf[p]:Tf[p]}
function rev(p,b){return b=="A"?Ar[p]:b=="C"?Cr[p]:b=="G"?Gr[p]:Tr[p]}
function cntb(p,b){return fwd(p,b)+rev(p,b)}
NR==FNR{Af[$1]=$2;Ar[$1]=$3;Cf[$1]=$4;Cr[$1]=$5;Gf[$1]=$6;Gr[$1]=$7;Tf[$1]=$8;Tr[$1]=$9;next}
{ gene=$1;ctg=$2;st=$3;canon=$4;plo=$5;pmid=$6;phi=$7;
  if(st=="+"){o[1]=plo;o[2]=pmid;o[3]=phi}else{o[1]=phi;o[2]=pmid;o[3]=plo}
  for(k=1;k<=3;k++){ p=ctg":"o[k]; dp=tot(p); if(dp<MINDP)continue;
    cb=substr(canon,k,1); eg=(st=="+")?cb:comp(cb);   # expected canonical genomic base
    split("A C G T",BB," ");
    for(z=1;z<=4;z++){ ab=BB[z]; if(ab==eg)continue; ac=cntb(p,ab); if(ac<MINALT)continue;
      mb=(st=="+")?ab:comp(ab);                        # alt mature base
      altac=canon; altac=substr(altac,1,k-1) mb substr(altac,k+1);
      vaf=100*ac/dp;
      printf "%s\t%s\t%s\t%s\t%s\t%d\t%s\t%.1f%%\t%d/%d\t%d\n", gene,ctg,st,canon,altac,o[k],ab"(mat "mb")",vaf,fwd(p,ab),rev(p,ab),dp;
    } } }
' pos_counts.tsv ac_table.bam.tsv | sort -t"$(printf '\t')" -k8,8nr > anticodon_variants.tsv

echo "===== ANTICODON-POSITION VARIANTS (genome-free; no VAF floor)  sample=$SAMPLE  minAltReads=$MINALT MAPQ>=$MINMAPQ baseQ>=20 ====="
printf "gene\tctg\tstr\tcanon\talt_anticodon\tpos\taltBase\tVAF\taltF/R\tdepth\n"
cat anticodon_variants.tsv
echo "---- $(wc -l < anticodon_variants.tsv) anticodon-position variants across $(wc -l < ac_table.bam.tsv) loci ----"

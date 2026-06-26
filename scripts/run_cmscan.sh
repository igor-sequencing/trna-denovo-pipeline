#!/usr/bin/env bash
set -euo pipefail
trap 'echo "CMSCAN EXIT $?"' EXIT
export PATH=/mnt/data/igor-shared/trnascan/env/bin:$PATH
CM=/mnt/data/igor-shared/trnascan/cmscan_SQ326CT3
cd "$CM"
THREADS=${THREADS:-24}
echo "[$(date)] cmscan START --rfam (threads=$THREADS, E<=1e-3)"
cmscan --rfam --cpu "$THREADS" -E 0.001 --noali \
  --tblout cmscan.tbl -o cmscan.out \
  trna_models.cm reads.fasta
echo "[$(date)] cmscan FINISHED"
echo "hit rows (E<=1e-3): $(grep -vc '^#' cmscan.tbl)"

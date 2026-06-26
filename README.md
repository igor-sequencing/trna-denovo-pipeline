# trna-denovo-pipeline

De novo assembly and anticodon analysis of human tRNAs from WGS data, using **NCBI SAUTE** (targeted de Bruijn assembler) + **tRNAscan-SE 2.0**, with a focus on screening for **non-canonical anticodons**.

The pipeline covers both **cytosolic** and **mitochondrial** tRNAs, starting from either raw FASTQ or an aligned BAM.

## Overview

```
                          ┌─ cytosolic: extract tRNA-locus reads ─┐
 BAM (or FASTQ) ──────────┤                                       ├─ SAUTE assembly ─ tRNAscan-SE ─ anticodons
                          └─ mito: extract ALL MT-contig reads ───┘                                  + non-canonical screen
```

1. **Reference set** (`build_ref.sh`) — cytosolic tRNAs from GtRNAdb hg38 (`hg38-tRNAs.fa`, 432 high-confidence; `hg38-tRNAs.bed`, 619 loci) + the 22 **mitochondrial** tRNAs built from the public **rCRS (NC_012920.1)** with `tRNAscan-SE -M mammal` (GtRNAdb deliberately excludes mito). Combined into `human_trna_ref.fa`. **No genome FASTA is needed** — only public tRNA databases + rCRS.
2. **Read retrieval** — two interchangeable front ends:
   - **FASTQ**: `retrieve_trna_reads.sh` — `bbduk.sh` k-mer match (k=25, both strands) against the tRNA reference. Streams the full WGS set in one pass.
   - **BAM** (faster, preferred when available): `samtools view -L <tRNA-loci.bed>` for cytosolic + `samtools view <MT>` for mitochondrial — index-driven, seconds.
3. **Assembly** — `saute` targeted to the tRNA reference.
4. **Validation / typing** — `tRNAscan-SE` on the assembled contigs (`-E` cytosolic, `-M mammal` mito) → isotype + anticodon calls.
5. **Non-canonical screen** (`scan_noncanonical_anticodons.sh`, genome-free) — reads each tRNA locus's 3 anticodon genomic positions (from the public `.ss` coordinates) directly from the BAM pileup and reports **every** anticodon-position alt allele supported by ≥3 quality reads (MAPQ≥20, baseQ≥20), with VAF / depth / strand. **No VAF cutoff** — real-vs-error is judged by read quality + strand balance, not allele frequency (a high-quality low-frequency allele is exactly the signal of interest). See Findings.

## Main entry point

```bash
./scripts/run_trna_pipeline.sh /path/to/sample.bam
```

Runs the whole thing for one BAM (cytosolic + mito + anticodon summary + non-canonical screen). Auto-detects the MT contig name (`MT`/`chrMT`/`chrM`/`M`) and reconciles GtRNAdb `chr`-prefixed coordinates to the BAM's contig naming.

## Hard-won parameters (read before changing)

- **SAUTE `--kmer 31 --secondary_kmer 21`** for cytosolic tRNAs. SAUTE's *automatic* kmer = ½ read length (75 for 150 bp reads), which **exceeds the length of most tRNAs** (reference min 59 bp; mito tRNAs 59–75 bp) — short tRNAs then produce no usable k-mers and silently fail to assemble. Forcing a small kmer fixes it (recovered 22/22 mito + full cytosolic isotype set vs only the long Leu/Ser/Tyr isotypes at kmer 75). SAUTE **requires both** `--kmer` and `--secondary_kmer` together.
- **Whole-mito assembly**: target the full 16.5 kb `MT` reference with **default (auto) kmer** — it's a long target, so the small-kmer rule does not apply.
- **`--max_variants 1000000`** uncaps per-tRNA variant output. ⚠️ Caveat: it also makes SAUTE emit alternative anticodon-position paths from read-level errors/paralogs, which can masquerade as "non-canonical" tRNAs. Always validate by pileup (see Findings).
- **cmscan has no `--fast`** — use `--rfam` for the fast preset (keeps the CM with aggressive filters); `--hmmonly` drops the CM entirely.
- **cmsearch/cmscan need FASTA** and are orders of magnitude slower than k-mer/alignment filtering — use them only as a confirmation step on the already-retrieved (tiny) read set, never on raw WGS.

## Scripts

| script | purpose |
|---|---|
| `build_ref.sh` | build combined cytosolic (GtRNAdb) + mito (`-M mammal`) tRNA reference |
| `run_trnascan_grch38.sh` | *optional* — scan a whole genome FASTA (`$GENOME`) from scratch; GtRNAdb already publishes this for GRCh38 |
| `run_trnascan.sh` | general parameterized tRNAscan-SE runner (`MODE` = `-E`, `-M mammal`, …; `$GENOME` input) |
| `retrieve_trna_reads.sh` | BBDuk k-mer retrieval of tRNA reads from FASTQ |
| `run_saute_k31.sh` | cytosolic SAUTE assembly (kmer 31) |
| `build_mt_saute.sh` | whole-mitochondrion: BAM `MT` extract → SAUTE (vs rCRS) → `tRNAscan-SE -M mammal` |
| `run_cmscan.sh` | Infernal cmscan confirmation of retrieved reads against tRNA CMs |
| `run_cyto_trnascan.sh` | `tRNAscan-SE -E` on cytosolic assembled contigs |
| `scan_noncanonical_anticodons.sh` | **genome-free** non-canonical anticodon screen from a BAM (pileup vs public `.ss` anticodon coords; no VAF floor) |
| `anticodon_stats.sh` | per-sample anticodon-count table (cytosolic + mito) from tRNAscan output |
| `run_trna_pipeline.sh` | **end-to-end** BAM → cytosolic + mito → anticodons + non-canonical screen |

## Environment

```bash
micromamba env create -f environment.yml      # or conda/mamba
```
tRNAscan-SE 2.0.12, Infernal 1.1.5, BBMap, samtools, seqkit, SKESA/SAUTE.

**No genome reference required** — the pipeline runs on the indexed **BAM** + the public **GtRNAdb** tRNA set/`.ss` + public **rCRS** + the **tRNAscan-SE** covariance models. (`run_trnascan_grch38.sh` is the only script that takes a whole-genome FASTA, via `$GENOME`, and it's optional.)

> Absolute paths inside the scripts (`/mnt/data/...`) point at the development environment (a Kubernetes pod with a persistent shared volume) and should be adjusted for other setups.

## Findings (GRCh38, 30× WGS)

- **Cytosolic**: ~49 distinct anticodons recovered from the assembled contigs, matching the 48 canonical GtRNAdb anticodons across all 20 amino acids + SeC + iMet.
- **Mitochondrial**: all **22/22** canonical mito tRNAs recovered from the whole-MT assembly.
- **Non-canonical screen**: the scanner reports every quality-supported anticodon-position alt allele with **no frequency cutoff** — real-vs-error is decided by base-Q / MAPQ / strand balance, not VAF, because a high-quality *low*-frequency allele (a minor gene copy, paralog, or mosaic) is exactly the signal of interest. Observed candidates are **A34 / A·NN** anticodons — the systematically "missing" class in Ehrlich et al. 2021, *Front. Mol. Biosci.* (PMC8007984): e.g. `Ala-AGC → ACC` at ~59% VAF (dominant; MAPQ 60, both strands) and `Asn-GTT → ATT` at ~9% VAF (MAPQ 60, baseQ 37, both strands — a genuine low-frequency candidate, **not** discarded as noise).

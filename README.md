# WES Germline Variant Calling Pipeline

A Nextflow DSL2 pipeline for whole exome sequencing (WES) germline variant discovery, joint genotyping, annotation and reporting across multi-family cohorts. Supports incremental batch processing via GenomicsDB UPDATE mode.

Author: Deepak Bharti, Clinical Bioinformatician, RCSI

contact: deepakbharti@rcsi.com

---

## Table of Contents

- [Overview](#overview)
- [Pipeline Architecture](#pipeline-architecture)
- [Tool Versions](#tool-versions)
- [Prerequisites](#prerequisites)
- [Sample Sheet Format](#sample-sheet-format)
- [Running the Pipeline](#running-the-pipeline)
  - [First Batch (CREATE mode)](#first-batch-create-mode)
  - [Second and Subsequent Batches (UPDATE mode)](#second-and-subsequent-batches-update-mode)
- [Output Structure](#output-structure)
- [Failure Recovery](#failure-recovery)
- [Key Parameters](#key-parameters)
- [Filter Thresholds](#filter-thresholds)

---

## Overview

```
FASTQ → QC → Alignment → BAM Processing → Variant Calling
                                               │
                              ┌────────────────┴────────────────┐
                         GATK HaplotypeCaller            DeepVariant
                              │                                  │
                         GenomicsDB                         GLnexus
                         (joint calling)              (joint calling)
                              │                                  │
                              └──────────┬──────────────────────┘
                                    Merge Callers
                                         │
                              Normalise (split multi-allelics)
                                         │
                              Variant Filtration (WES hard filters)
                                         │
                                   SpliceAI Annotation
                                         │
                                  Per-sample ANNOVAR
                                         │
                                   Excel Report
```

---

## Pipeline Architecture

| Stage | Process | Description |
|---|---|---|
| QC | FastQC, Fastp | Raw read quality control and trimming |
| Alignment | BWA-MEM2 | Read alignment to hg38 reference |
| BAM Processing | MarkDuplicates, BQSR | Duplicate marking, base quality recalibration |
| Variant Calling | HaplotypeCaller (ploidy-aware) | Per-sample GVCF generation |
| Variant Calling | DeepVariant | CNN-based per-sample variant calling |
| Joint Calling | GenomicsDB + GenotypeGVCFs | Incremental joint genotyping (GATK) |
| Joint Calling | GLnexus | Joint genotyping (DeepVariant) |
| Merge | BCFtools | Caller concordance tagging (BOTH/GATK_ONLY/DV_ONLY) |
| Normalisation | BCFtools norm | Split multi-allelics, left-align indels |
| Filtering | GATK VariantFiltration | WES-tuned hard filters (SNP + INDEL) |
| Annotation | SpliceAI | Splice site impact scoring (deep learning) |
| Annotation | ANNOVAR | Functional annotation (refGene, gnomAD, ClinVar etc.) |
| Reporting | Python | Per-sample Excel workbooks |
| QC | MultiQC | Aggregated QC report |

---

## Tool Versions

| Tool | Version | Purpose |
|---|---|---|
| Nextflow | ≥ 25.04 | Pipeline orchestration |
| BWA-MEM2 | 2.2.1 | Read alignment |
| SAMtools | 1.22 | BAM processing |
| GATK | 4.6.2.0 | HaplotypeCaller, GenomicsDB, GenotypeGVCFs, VariantFiltration |
| GenomicsDB | 1.5.5 | Incremental joint calling database |
| DeepVariant | 1.6.1 | CNN variant calling |
| GLnexus | 1.4.1 | DeepVariant joint calling |
| BCFtools | 1.23 | VCF manipulation, merging, normalisation |
| Fastp | 0.23.4 | Read trimming and QC |
| FastQC | 0.12.1 | Read quality assessment |
| Mosdepth | 0.3.9 | Coverage analysis |
| SpliceAI | 1.3.1 | Splice site variant scoring |
| ANNOVAR | 2020-06-08 | Functional variant annotation |
| Peddy | 0.4.8 | Pedigree QC and sex/ancestry inference |
| MultiQC | 1.25 | Aggregated QC reporting |
| Python | 3.10+ | Post-processing and Excel report generation |
| Singularity | ≥ 3.8 | Container runtime |

---

## Prerequisites

- SLURM HPC cluster with Singularity ≥ 3.8
- Nextflow ≥ 25.04
- Reference files: hg38 FASTA, dbSNP, known indels, ANNOVAR databases
- Singularity image directory configured in `nextflow.config`

---

## Sample Sheet Format

The pipeline takes a CSV sample sheet. Each row represents one sample member.

### Column Definitions

| Column | Description | Values |
|---|---|---|
| `sample_id` | Unique sample identifier | e.g. `FAM01_PROBAND` |
| `family_id` | Family identifier | e.g. `FAM01` |
| `fastq_r1` | Full path to R1 FASTQ | `/data/FAM01_R1.fastq.gz` |
| `fastq_r2` | Full path to R2 FASTQ | `/data/FAM01_R2.fastq.gz` |
| `sex` | Biological sex | `1` or `2` |
| `phenotype` | phenotype status | `2` (phenotype) or `1` or `0` or `-9`(unphenotype) |
| `paternal_id` | Father's sample ID | `FAM01_FATHER` or `0` if unknown |
| `maternal_id` | Mother's sample ID | `FAM01_MOTHER` or `0` if unknown |

* 1 = Male, 2 = Female
### Example `sample.csv`

```csv
sample_id,family_id,fastq_r1,fastq_r2,sex,phenotype,paternal_id,maternal_id
FAM01_PROBAND,FAM01,/data/fastq/FAM01_PRB_R1.fastq.gz,/data/fastq/FAM01_PRB_R2.fastq.gz,1,2,FAM01_FATHER,FAM01_MOTHER
FAM01_FATHER,FAM01,/data/fastq/FAM01_FAT_R1.fastq.gz,/data/fastq/FAM01_FAT_R2.fastq.gz,2,2,0,0
FAM01_MOTHER,FAM01,/data/fastq/FAM01_MOT_R1.fastq.gz,/data/fastq/FAM01_MOT_R2.fastq.gz,1,1,0,0
FAM02_PROBAND,FAM02,/data/fastq/FAM02_PRB_R1.fastq.gz,/data/fastq/FAM02_PRB_R2.fastq.gz,2,2,FAM02_FATHER,FAM02_MOTHER
FAM02_FATHER,FAM02,/data/fastq/FAM02_FAT_R1.fastq.gz,/data/fastq/FAM02_FAT_R2.fastq.gz,1,1,0,0
FAM02_MOTHER,FAM02,/data/fastq/FAM02_MOT_R1.fastq.gz,/data/fastq/FAM02_MOT_R2.fastq.gz,2,1,0,0
```

> **Note:** Trio structure (proband + both parents) is required for de novo variant calling. Duos and singletons are supported but de novo calling will be skipped.

---

## Running the Pipeline

The pipeline is submitted via the **`run_pipeline.slurm`** master script. All batch configuration is done by editing variables at the **top of this script** - no separate parameter file is required.

### The Master Submit Script

The top of `run_pipeline.slurm` contains all user-configurable settings:

```bash
# ============================================================================
# USER CONFIGURATION - EDIT THESE FOR EACH BATCH
# ============================================================================

BATCH_NAME="COHORT_BATCH1"        # Batch identifier - used in all output filenames
INPUT_CSV="batch1_samples.csv"    # Sample sheet for THIS batch only
```

The script automatically passes all settings to Nextflow and handles:
- GenomicsDB CREATE vs UPDATE mode
- Automatic backup of GenomicsDB before UPDATE runs
- Master PED file management across batches
- Correct output namespacing per batch

---

### First Batch (CREATE mode)

**1. Prepare sample sheet** - see [Sample Sheet Format](#sample-sheet-format)

**2. Edit the top of `run_pipeline.slurm`:**

```bash
BATCH_NAME="COHORT_BATCH1"
INPUT_CSV="batch1_samples.csv"
```

**3. Submit:**

```bash
sbatch run_pipeline.slurm
```

**4. Monitor:**

```bash
# Watch jobs
watch squeue -u $USER

# Watch Nextflow log live
tail -f .nextflow.log | grep -E "Submitted|COMPLETED|FAILED|ERROR"
```

---

### Second and Subsequent Batches (UPDATE mode)

Each subsequent batch adds new samples to the existing GenomicsDB and re-genotypes the **entire cohort jointly**. The joint VCF for BATCH2 contains all BATCH1 + BATCH2 samples combined.

**1. Prepare sample sheet containing only the NEW batch samples.** Do not re-include previous batch samples - the pipeline reads existing samples from GenomicsDB automatically.

**2. Edit the top of `run_pipeline.slurm`:**

```bash
BATCH_NAME="COHORT_BATCH2"
INPUT_CSV="batch2_samples.csv"
```

**3. Submit:**

```bash
sbatch run_pipeline.slurm
```

> **Backup:** The script automatically backs up the existing GenomicsDB before the UPDATE runs. The backup path is printed in the SLURM output log. See [Failure Recovery](#failure-recovery) for restore instructions.

---

## Output Structure

All outputs are namespaced by `BATCH_NAME` - results from different batches never overwrite each other.

```
cohort_results/
├── variants/
│   ├── gvcf/                              # Per-sample GVCFs (HaplotypeCaller)
│   ├── joint/
│   │   └── COHORT_BATCH1/
│   │       ├── COHORT_BATCH1.joint.vcf.gz          # Joint genotyped (all samples)
│   │       └── COHORT_BATCH1.joint.vcf.gz.tbi
│   ├── merged/
│   │   └── COHORT_BATCH1/
│   │       ├── COHORT_BATCH1.merged.caller_tagged.vcf.gz   # GATK + DV merged
│   │       └── COHORT_BATCH1.concordance_stats.txt
│   ├── gatk_filtered/
│   │   └── COHORT_BATCH1/
│   │       ├── COHORT_BATCH1.gatk_filtered.vcf.gz   # GATK-only hard filtered
│   │       └── COHORT_BATCH1.gatk_filtered.vcf.gz.tbi
│   ├── spliceai/
│   │   └── COHORT_BATCH1/
│   │       ├── COHORT_BATCH1.filtered.spliceai.merged.vcf.gz
│   │       └── COHORT_BATCH1.spliceai_summary.txt
│   ├── deepvariant/                       # Per-sample DeepVariant VCFs
│   └── glnexus/                           # GLnexus joint VCF
├── genomicsdb/
│   └── genomicsdb/                        # GenomicsDB workspace (used for UPDATE)
├── annotation/
│   └── <SAMPLE_ID>/
│       └── <SAMPLE_ID>.annotated.xlsx     # Per-sample Excel report
├── qc/
│   ├── fastqc/
│   ├── fastp/
│   ├── coverage/
│   └── multiqc_report.html
└── pipeline_info/
    ├── report_<timestamp>.html
    ├── timeline_<timestamp>.html
    └── trace_<timestamp>.txt
```

---

## Failure Recovery

If the pipeline fails or is cancelled (e.g. via `scancel`), follow this **exact sequence** before resuming:

### Step 1 - Identify the failure

```bash
grep -E "ERROR|FAILED|Cause" .nextflow.log | tail -20
```

### Step 2 - Remove the master PED file

The pipeline auto-generates a merged pedigree (`master.ped`) across batches. If interrupted mid-run, a partial PED file will cause failures on resume:

```bash
# Remove the master PED - it is regenerated automatically on resume
rm -f /path/to/results/ped/master.ped
```

### Step 3 - Remove the Nextflow LOCK file

The LOCK file prevents two Nextflow runs from corrupting the cache. It is not released cleanly after `scancel`:

```bash
rm -f /path/to/workdir/.nextflow/cache/*/db/LOCK
```

### Step 4 - Check for runaway files

Occasionally a failed process can write a very large partial output file:

```bash
# Check for files > 10GB in the work directory
find /path/to/workdir -size +10G -newer .nextflow.log
# Remove any unexpected large files before resuming
```

### Step 5 - Resume

Add `-resume` to the Nextflow command in `run_pipeline.slurm` then resubmit:

```bash
sbatch run_pipeline.slurm
```

> `-resume` uses Nextflow's content-hash cache - only failed or new tasks rerun. All successfully completed tasks are skipped instantly.

---

### Restoring GenomicsDB from Backup (UPDATE runs only)

If an UPDATE run fails and corrupts the GenomicsDB, restore from the automatic backup:

```bash
# The backup path is printed in the SLURM output log at the start of each UPDATE run
# Example: "Backup created at: /path/to/results/backups/COHORT_BATCH2_20260330_113000"

# Restore
rm -rf /path/to/results/genomicsdb/genomicsdb
cp -rL /path/to/results/backups/COHORT_BATCH2_<timestamp>/genomicsdb \
        /path/to/results/genomicsdb/genomicsdb
```

---

## Key Parameters

All parameters are set in `run_pipeline.slurm`. The following are edited per batch:

| SLURM variable | Nextflow param | Description |
|---|---|---|
| `BATCH_NAME` | `--batch_name` | Batch identifier - used in all output filenames |
| `INPUT_CSV` | `--input` | Sample CSV for this batch (new samples only for UPDATE) |
| `RUN_MODE` | `--run_mode` | `CREATE` (first batch) or `UPDATE` (subsequent) |
| `OUTDIR` | `--outdir` | Results output directory |

---

## Filter Thresholds

### SNP Hard Filters (GATK VariantFiltration, WES-tuned)

| Filter name | Expression | Rationale |
|---|---|---|
| `QD_filter` | `QD < 2.0` | Low quality normalised by depth |
| `FS_filter` | `FS > 60.0` | Strand bias (Fisher's Exact Test) |
| `MQ_filter` | `MQ < 40.0` | Poor mapping quality |
| `SOR_filter` | `SOR > 3.0` | Strand odds ratio bias |
| `MQRankSum_filter` | `MQRankSum < -12.5` | Mapping quality rank sum |
| `ReadPosRankSum_filter` | `ReadPosRankSum < -8.0` | Read position rank sum |
| `QUAL_filter` | `QUAL < 30.0` | Raw call confidence |
| `SnpCluster` | 3 SNPs within 25bp | Clustered SNP artifact - 25bp window WES-tuned |

### INDEL Hard Filters

| Filter name | Expression | Rationale |
|---|---|---|
| `QD_filter` | `QD < 2.0` | Low quality normalised by depth |
| `FS_filter` | `FS > 200.0` | Strand bias (relaxed for INDELs) |
| `SOR_filter` | `SOR > 10.0` | Strand odds ratio (relaxed for INDELs) |
| `QUAL_filter` | `QUAL < 30.0` | Raw call confidence |

> **Note:** `ReadPosRankSum` and `MQ` filters are intentionally omitted for INDELs per GATK best practice. INDEL alignment mechanics produce naturally worse position bias scores unrelated to variant quality.

---

## Notes

- **Dual caller concordance:** Variants called by both GATK and DeepVariant (`CALLER=BOTH`) have the highest confidence. `CALLER=GATK_ONLY` or `CALLER=DV_ONLY` variants warrant additional scrutiny.
- **Multi-allelic normalisation:** Sites with multiple ALT alleles are split into biallelic records before filtering to prevent filter expressions being silently skipped.
- **PASS-only annotation:** Only `FILTER=PASS` variants are passed to per-sample annotation. Hard-filtered variants are excluded from Excel reports.
- **SpliceAI efficiency:** SpliceAI runs on the cohort-level VCF before per-sample splitting - scores are computed once per unique variant position regardless of how many samples carry it.
- **Batch isolation:** Each batch's outputs are namespaced by `BATCH_NAME`. Running BATCH2 never overwrites BATCH1 results.
=======

# WES/WGS Germline Variant Calling Pipeline v3.0

## Quick Start

### 1. Deploy Pipeline

```bash
# Copy pipeline to your project directory
cp -r wes-pipeline /path/to/wes_pipeline_v2
cd /path/to/wes_pipeline_v2
```

### 2. Prepare Sample CSV

Create your sample sheet (see `sample_sheet_template.csv` for format):

```csv
sample_id,family_id,fastq_r1,fastq_r2,sex,phenotype,paternal_id,maternal_id
FAM01_PROBAND,FAM01,/data/fastq/FAM01_PRB_R1.fastq.gz,/data/fastq/FAM01_PRB_R2.fastq.gz,1,2,FAM01_FATHER,FAM01_MOTHER
FAM01_FATHER,FAM01,/data/fastq/FAM01_FAT_R1.fastq.gz,/data/fastq/FAM01_FAT_R2.fastq.gz,2,2,0,0
FAM01_MOTHER,FAM01,/data/fastq/FAM01_MOT_R1.fastq.gz,/data/fastq/FAM01_MOT_R2.fastq.gz,1,1,0,0
FAM02_PROBAND,FAM02,/data/fastq/FAM02_PRB_R1.fastq.gz,/data/fastq/FAM02_PRB_R2.fastq.gz,2,2,FAM02_FATHER,FAM02_MOTHER
FAM02_FATHER,FAM02,/data/fastq/FAM02_FAT_R1.fastq.gz,/data/fastq/FAM02_FAT_R2.fastq.gz,1,1,0,0
FAM02_MOTHER,FAM02,/data/fastq/FAM02_MOT_R1.fastq.gz,/data/fastq/FAM02_MOT_R2.fastq.gz,2,1,0,0
```

### 3. Edit SLURM Script

Edit `run_pipeline.slurm`:

```bash
# Line 24-28: Edit these for each batch
BATCH_NAME="batch1"
INPUT_CSV="sample_sheet_batch1.csv"
```

### 4. Run Pipeline

```bash
# First batch - creates GenomicsDB and master.ped
sbatch run_pipeline.slurm

# Check status
squeue -u $USER

# Monitor logs
tail -f nextflow_*.out
```

### 5. Subsequent Batches

```bash
# Edit script for new batch
vi run_pipeline.slurm
# Change: BATCH_NAME="batch2"
# Change: INPUT_CSV="sample_sheet_batch2.csv"

# Run - auto-detects incremental mode
sbatch run_pipeline.slurm
```

---

## Output Structure

```
results_v2/
├── annotation/per_sample/
│                    ├── sm1.hg38_multianno.csv    # ← Main output with GT, Caller, DeNovo columns
│                    └── sm1.hg38.avinput        
├── variants/
│   ├── joint/
│   │   └── joint.vcf.gz             # Joint-called VCF (all samples)
│   ├── gvcf/                        # Per-sample GVCFs
│   └── merged/                      # Merged GATK+DV VCF (if DeepVariant enabled)
├── bam/
│   ├── markdup/                     # MarkDuplicates BAMs
│   └── bqsr/                        # BQSR-recalibrated BAMs
├── ped/
│   ├── batch.ped                    # Current batch PED
│   └── master.ped                   # Cumulative PED (all samples)
├── qc/
│   ├── fastqc/
│   ├── fastp/
│   ├── coverage/
│   ├── peddy/
│   └── multiqc/
│       └── multiqc_report.html      # Aggregated QC report
├── genomicsdb/                      # GenomicsDB (all samples)
└── pipeline_info/
    ├── timeline.html
    ├── report.html
    └── trace.txt
```

---

## Key Output: multianno.csv Columns

| Column | Description |
|--------|-------------|
| Chr, Start, End, Ref, Alt | Variant coordinates |
| Func.refGene | Gene function |
| Gene.refGene | Gene name |
| ExonicFunc.refGene | Exonic function (missense, nonsense, etc.) |
| AAChange.refGene | Amino acid change |
| gnomad41_exome | gnomAD exome allele frequency |
| clinvar_20250721 | ClinVar annotation |
| **Caller** | `BOTH`, `GATK_ONLY`, `DV_ONLY`, or `GATK` |
| **DeNovo** | `hiConfDeNovo:sample`, `loConfDeNovo:sample`, or `.` |
| **GT_sample1** | `het`, `hom-alt`, `hom-ref`, or `missing` |
| **GT_sample2** | (one column per batch sample) |

---

## Enable DeepVariant (Optional)

```bash
# In run_pipeline.slurm, change:
RUN_DEEPVARIANT="true"
```

This will:
1. Run DeepVariant on MarkDuplicates BAMs (not BQSR)
2. Run GLnexus for DeepVariant joint calling
3. Merge GATK + DeepVariant VCFs
4. Tag each variant: `CALLER=BOTH/GATK_ONLY/DV_ONLY`

---

## Troubleshooting

### Check Logs
```bash
# Nextflow log
cat .nextflow.log

# SLURM output
cat nextflow_<jobid>.out

# Failed task work directory
cat work/<hash>/<hash>/.command.err
```

### Resume Failed Run
```bash
# Pipeline automatically resumes with -resume flag
sbatch run_pipeline.slurm
```

### Force Fresh Start
```bash
# Remove -resume from script, or:
rm -rf work/
sbatch run_pipeline.slurm
```

---

## Container Verification

Ensure these containers exist in `/path/to/docker_iamges`:

```
fastqc_0.12.1.sif
fastp_0.23.4.sif
bwa-mem2-samtool2.3.1_1.22.1.sif
samtools_1.23.sif
gatk4_4.6.2.0.sif
deepvariant_1.9.0.sif
glnexus_1.2.7.sif          # ← Required for DeepVariant mode
bcftool_1.23.sif
mosdepth_0.3.11.sif
multiqc_1.33.sif
peddy_0.4.8.sif
annovar_2020Jun08.sif
python_3.14.2.sif
```
=======
# WES_germaline_variant_calling
Containerized Nextflow pipeline for WES, easily extendable to WGS (via VQSR). It integrates FastQC/MultiQC for QC, fastp for preprocessing, BWA-MEM2 for alignment, GATK 4.2.6 for variant calling, and ANNOVAR for annotation, with optional DeepVariant and SpliceAI support.


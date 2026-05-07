# WES Germline Variant Calling Pipeline v3.0

A Nextflow DSL2 pipeline for whole exome sequencing (WES) germline variant discovery, joint genotyping, annotation and reporting across multi-family cohorts. Supports incremental batch processing via GenomicsDB UPDATE mode.

**Author:** Deepak Bharti, Clinical Bioinformatician, RCSI  
**Contact:** deepakbharti@rcsi.com

---

## Table of Contents

- [Overview](#overview)
- [Key Improvements in v3.0](#key-improvements-in-v30)
- [Pipeline Architecture](#pipeline-architecture)
- [Tool Versions](#tool-versions)
- [Prerequisites](#prerequisites)
- [Sample Sheet Format](#sample-sheet-format)
- [Running the Pipeline](#running-the-pipeline)
- [Output Structure](#output-structure)
- [Failure Recovery](#failure-recovery)
- [Key Parameters](#key-parameters)
- [Filter Thresholds](#filter-thresholds)

---

## Overview

```
FASTQ (single or multi-lane) → QC → Alignment (per lane) → BAM Merge (if multi-lane)
                                                                      │
                                                              BAM Processing
                                                         (MarkDuplicates + BQSR)
                                                                      │
                                         ┌────────────────────────────┴────────────────────────────┐
                                  GATK HaplotypeCaller                                      DeepVariant
                                         │                                                          │
                                    GenomicsDB                                                 GLnexus
                                  (joint calling)                                          (joint calling)
                                         │                                                          │
                                         └──────────────────────┬─────────────────────────────────┘
                                                           Merge Callers
                                                                 │
                                                    DP=0 Masking (BCFTOOLS SETGT)
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

## Key Improvements in v3.0

### 1. Lane-aware FASTQ Support
Multiple FASTQ pairs per sample (e.g. L001, L002 from the same flowcell) are now natively supported. Each lane is aligned independently with a lane-specific read group tag, then merged into a single per-sample BAM before downstream processing.

- `lane` column added to samplesheet (optional  defaults to `L001` when absent)
- Per-lane output naming prevents filename collisions (`sample.L001.sorted.bam`)
- Read group tags include `ID=sample.lane`, `PU=lane`, `SM=sample` for correct BQSR grouping
- Single-lane samples pass through with zero overhead  no merge step is invoked

### 2. PED File Deduplication
Multi-lane samplesheets (multiple rows per sample) no longer produce duplicate PED entries. Each sample produces exactly one PED line regardless of lane count, preventing downstream peddy and de novo calling errors.

### 3. Per-lane QC Naming
FastQC and Fastp outputs are named `sample.lane.*` to prevent MultiQC staging collisions when the same sample has multiple lanes. All per-lane reports are aggregated correctly in the MultiQC report.

### 4. DP=0 Genotype Masking (Multi-kit Cohort Support)
A new `BCFTOOLS_SETGT` step runs after joint genotyping and converts `DP=0` genotypes from `0/0` to `./.` (missing). This prevents false de novo calls and allele frequency distortion in cohorts where samples were captured with different exome kits (preferably from same provider). Samples with no reads in kit-exclusive regions are correctly treated as missing rather than homozygous reference. This step has no effect when all samples share the same capture kit.

### 5. Merged BAM Persistent Storage
Merged BAMs (multi-lane samples only) are published to a configurable persistent path (`params.merged_bam_path`) for downstream reuse without reprocessing.

---

## Pipeline Architecture

| Stage | Process | Description |
|---|---|---|
| QC | FastQC, Fastp | Per-lane read quality control and trimming |
| Alignment | BWA-MEM2 | Per-lane alignment to hg38 with lane-aware RG tags |
| BAM Merge | SAMtools merge | Merge per-lane BAMs into per-sample BAM (multi-lane only) |
| BAM Processing | MarkDuplicates, BQSR | Duplicate marking, base quality recalibration |
| Variant Calling | HaplotypeCaller (ploidy-aware) | Per-sample GVCF generation |
| Variant Calling | DeepVariant | CNN-based per-sample variant calling |
| Joint Calling | GenomicsDB + GenotypeGVCFs | Incremental joint genotyping (GATK) |
| Joint Calling | GLnexus | Joint genotyping (DeepVariant) |
| Masking | BCFtools setGT | DP=0 → ./. for multi-kit cohort safety |
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
| SAMtools | 1.22 | BAM processing and merging |
| GATK | 4.6.2.0 | HaplotypeCaller, GenomicsDB, GenotypeGVCFs, VariantFiltration |
| GenomicsDB | 1.5.5 | Incremental joint calling database |
| DeepVariant | 1.9.0 | CNN variant calling |
| GLnexus | 1.2.7 | DeepVariant joint calling |
| BCFtools | 1.23 | VCF manipulation, merging, normalisation, setGT masking |
| Fastp | 0.23.4 | Read trimming and QC |
| FastQC | 0.12.1 | Read quality assessment |
| Mosdepth | 0.3.11 | Coverage analysis |
| SpliceAI | 1.3.1 | Splice site variant scoring |
| ANNOVAR | 2020-06-08 | Functional variant annotation |
| Peddy | 0.4.8 | Pedigree QC and sex/ancestry inference |
| MultiQC | 1.33 | Aggregated QC reporting |
| Python | 3.14.2 | Post-processing and Excel report generation |
| Singularity | ≥ 3.8 | Container runtime |

---

## Prerequisites

- SLURM HPC cluster with Singularity ≥ 3.8
- Nextflow ≥ 25.04
- Reference files: hg38 FASTA, dbSNP, known indels, ANNOVAR databases
- Singularity image directory configured in `nextflow.config`

---

## Sample Sheet Format

The pipeline takes a CSV sample sheet. Each row represents one sample (or one lane of a sample).

### Column Definitions

| Column | Required | Description | Values |
|---|---|---|---|
| `family_id` | Yes | Family identifier | e.g. `FAM01` |
| `sample_id` | Yes | Unique sample identifier | e.g. `FAM01_PROBAND` |
| `lane` | No | Lane identifier  omit for single-lane | `L001`, `L002` etc. Defaults to `L001` |
| `paternal_id` | No | Father's sample ID | `FAM01_FATHER` or `0` |
| `maternal_id` | No | Mother's sample ID | `FAM01_MOTHER` or `0` |
| `sex` | No | Biological sex | `1`=Male, `2`=Female, `0`=Unknown |
| `phenotype` | No | Phenotype status | `2`=Affected, `1`=Unaffected, `0`/`-9`=Unknown |
| `fastq_r1` | Yes | Full path to R1 FASTQ | `/data/FAM01_R1.fastq.gz` |
| `fastq_r2` | Yes | Full path to R2 FASTQ | `/data/FAM01_R2.fastq.gz` |

> **Column order does not matter**  the pipeline reads by column name, not position.

> **Note:** Trio structure (proband + both parents) is required for de novo variant calling. Duos and singletons are supported but de novo calling will be skipped. If sex information is unavailable, `0`/`-9` is accepted but X/Y chromosome calls are unreliable.

---

### Example: Standard Samplesheet (No Lane Column)

Use this format when each sample has a single pair of FASTQ files.

```csv
family_id,sample_id,paternal_id,maternal_id,sex,phenotype,fastq_r1,fastq_r2
FAM01,FAM01_PROBAND,FAM01_FATHER,FAM01_MOTHER,1,2,/data/fastq/FAM01_PRB_R1.fastq.gz,/data/fastq/FAM01_PRB_R2.fastq.gz
FAM01,FAM01_FATHER,0,0,1,1,/data/fastq/FAM01_FAT_R1.fastq.gz,/data/fastq/FAM01_FAT_R2.fastq.gz
FAM01,FAM01_MOTHER,0,0,2,1,/data/fastq/FAM01_MOT_R1.fastq.gz,/data/fastq/FAM01_MOT_R2.fastq.gz
FAM02,FAM02_PROBAND,FAM02_FATHER,FAM02_MOTHER,2,2,/data/fastq/FAM02_PRB_R1.fastq.gz,/data/fastq/FAM02_PRB_R2.fastq.gz
FAM02,FAM02_FATHER,0,0,1,1,/data/fastq/FAM02_FAT_R1.fastq.gz,/data/fastq/FAM02_FAT_R2.fastq.gz
FAM02,FAM02_MOTHER,0,0,2,1,/data/fastq/FAM02_MOT_R1.fastq.gz,/data/fastq/FAM02_MOT_R2.fastq.gz
```

Each sample has one row → aligned once → no BAM merge step.

---

### Example: Lane-aware Samplesheet (Multi-lane)

Use this format when a sample was sequenced across multiple lanes (e.g. high-depth sequencing split across L001 and L002).

```csv
family_id,sample_id,lane,paternal_id,maternal_id,sex,phenotype,fastq_r1,fastq_r2
FAM01,FAM01_PROBAND,L001,FAM01_FATHER,FAM01_MOTHER,1,2,/data/fastq/FAM01_PRB_L001_R1.fastq.gz,/data/fastq/FAM01_PRB_L001_R2.fastq.gz
FAM01,FAM01_PROBAND,L002,FAM01_FATHER,FAM01_MOTHER,1,2,/data/fastq/FAM01_PRB_L002_R1.fastq.gz,/data/fastq/FAM01_PRB_L002_R2.fastq.gz
FAM01,FAM01_FATHER,L001,0,0,1,1,/data/fastq/FAM01_FAT_L001_R1.fastq.gz,/data/fastq/FAM01_FAT_L001_R2.fastq.gz
FAM01,FAM01_MOTHER,L001,0,0,2,1,/data/fastq/FAM01_MOT_L001_R1.fastq.gz,/data/fastq/FAM01_MOT_L001_R2.fastq.gz
```

- `FAM01_PROBAND` has two rows (L001 + L002) → aligned separately → BAMs merged → single merged BAM enters BAM processing
- `FAM01_FATHER` and `FAM01_MOTHER` have one row each → pass through directly, no merge
- PED file contains exactly one entry per `sample_id` regardless of lane count
- Mixing single-lane and multi-lane samples in the same samplesheet is fully supported

---

### Example: Quad Family (Two Affected Siblings, Shared Parents)

```csv
family_id,sample_id,paternal_id,maternal_id,sex,phenotype,fastq_r1,fastq_r2
FAM03,FAM03_FATHER,0,0,1,1,/data/FAM03_DAD_R1.fastq.gz,/data/FAM03_DAD_R2.fastq.gz
FAM03,FAM03_MOTHER,0,0,2,1,/data/FAM03_MOM_R1.fastq.gz,/data/FAM03_MOM_R2.fastq.gz
FAM03,FAM03_CHILD1,FAM03_FATHER,FAM03_MOTHER,1,2,/data/FAM03_CH1_R1.fastq.gz,/data/FAM03_CH1_R2.fastq.gz
FAM03,FAM03_CHILD2,FAM03_FATHER,FAM03_MOTHER,1,2,/data/FAM03_CH2_R1.fastq.gz,/data/FAM03_CH2_R2.fastq.gz
```

Both children share the same parents  valid standard PED quad structure. Each child is evaluated independently against the same parental genotypes for de novo calling.

---

## Running the Pipeline

The pipeline is submitted via the **`run_pipeline.slurm`** master script. All batch configuration is done by editing variables at the top of this script.

### First Batch (CREATE mode)

**1. Prepare sample sheet**  see [Sample Sheet Format](#sample-sheet-format)

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
watch squeue -u $USER
tail -f .nextflow.log | grep -E "Submitted|COMPLETED|FAILED|ERROR"
```

---

### Second and Subsequent Batches (UPDATE mode)

Each subsequent batch adds new samples to the existing GenomicsDB and re-genotypes the entire cohort jointly. The joint VCF for BATCH2 contains all BATCH1 + BATCH2 samples combined.

**1. Prepare sample sheet containing only the NEW batch samples.**

**2. Edit the top of `run_pipeline.slurm`:**

```bash
BATCH_NAME="COHORT_BATCH2"
INPUT_CSV="batch2_samples.csv"
```

**3. Submit:**

```bash
sbatch run_pipeline.slurm
```

> **Backup:** The script automatically backs up the existing GenomicsDB before the UPDATE runs.

---

## Output Structure

```
cohort_results/
├── bam/
│   └── merged/                            # Per-sample merged BAMs (multi-lane samples only)
│       ├── SAMPLE001.merged.bam
│       └── SAMPLE001.merged.bam.bai
├── variants/
│   ├── gvcf/                              # Per-sample GVCFs (HaplotypeCaller)
│   ├── joint/
│   │   └── COHORT_BATCH1/
│   │       ├── COHORT_BATCH1.joint.vcf.gz
│   │       └── COHORT_BATCH1.joint.vcf.gz.tbi
│   ├── merged/
│   │   └── COHORT_BATCH1/
│   │       ├── COHORT_BATCH1.merged.caller_tagged.vcf.gz
│   │       └── COHORT_BATCH1.concordance_stats.txt
│   ├── gatk_filtered/
│   │   └── COHORT_BATCH1/
│   │       ├── COHORT_BATCH1.gatk_filtered.vcf.gz
│   │       └── COHORT_BATCH1.gatk_filtered.vcf.gz.tbi
│   ├── spliceai/
│   │   └── COHORT_BATCH1/
│   │       ├── COHORT_BATCH1.filtered.spliceai.merged.vcf.gz
│   │       └── COHORT_BATCH1.spliceai_summary.txt
│   ├── deepvariant/                       # Per-sample DeepVariant GVCFs
│   └── glnexus/                           # GLnexus joint VCF
├── genomicsdb/
│   └── genomicsdb/                        # GenomicsDB workspace
├── ped/
│   ├── batch.ped                          # Current batch PED (deduplicated)
│   ├── master.ped                         # Cumulative PED (all batches)
│   └── master.ped_backup_<timestamp>      # Auto-backup before each UPDATE
├── annotation/
│   └── <SAMPLE_ID>/
│       └── <SAMPLE_ID>.annotated.xlsx     # Per-sample Excel report
├── qc/
│   ├── fastqc/                            # Per-lane FastQC reports
│   ├── fastp/                             # Per-lane Fastp reports
│   ├── coverage/
│   ├── peddy/
│   └── multiqc_report.html
└── pipeline_info/
    ├── report_<timestamp>.html
    ├── timeline_<timestamp>.html
    └── trace_<timestamp>.txt
```

---

## Failure Recovery

If the pipeline fails or is cancelled, follow this exact sequence before resuming:

### Step 1  Identify the failure

```bash
grep -E "ERROR|FAILED|Cause" .nextflow.log | tail -20
```

### Step 2  Remove the master PED file

```bash
rm -f /path/to/results/ped/master.ped
```

### Step 3  Remove the Nextflow LOCK file

```bash
rm -f /path/to/workdir/.nextflow/cache/*/db/LOCK
```

### Step 4  Check for runaway files

```bash
find /path/to/workdir -size +10G -newer .nextflow.log
```

### Step 5  Resume

```bash
sbatch run_pipeline.slurm   # -resume is set automatically
```

### Restoring GenomicsDB from Backup (UPDATE runs only)

```bash
rm -rf /path/to/results/genomicsdb/genomicsdb
cp -rL /path/to/results/backups/COHORT_BATCH2_<timestamp>/genomicsdb \
        /path/to/results/genomicsdb/genomicsdb
```

---

## Key Parameters

| SLURM variable | Nextflow param | Description |
|---|---|---|
| `BATCH_NAME` | `--batch_name` | Batch identifier used in all output filenames |
| `INPUT_CSV` | `--input` | Sample CSV for this batch |
| `RUN_MODE` | `--run_mode` | `CREATE` (first batch) or `UPDATE` (subsequent) |
| `OUTDIR` | `--outdir` | Results output directory |
| `merged_bam_path` | `--merged_bam_path` | Persistent path for merged BAMs (multi-lane) |

---

## Filter Thresholds

### SNP Hard Filters (WES-tuned)

| Filter | Expression | Rationale |
|---|---|---|
| `QD_filter` | `QD < 2.0` | Low quality normalised by depth |
| `FS_filter` | `FS > 60.0` | Strand bias |
| `MQ_filter` | `MQ < 40.0` | Poor mapping quality |
| `SOR_filter` | `SOR > 3.0` | Strand odds ratio |
| `MQRankSum_filter` | `MQRankSum < -12.5` | Mapping quality rank sum |
| `ReadPosRankSum_filter` | `ReadPosRankSum < -8.0` | Read position rank sum |
| `QUAL_filter` | `QUAL < 30.0` | Raw call confidence |
| `SnpCluster` | 3 SNPs within 35bp | Clustered SNP artifact |

### INDEL Hard Filters

| Filter | Expression | Rationale |
|---|---|---|
| `QD_filter` | `QD < 2.0` | Low quality normalised by depth |
| `FS_filter` | `FS > 200.0` | Strand bias (relaxed for INDELs) |
| `SOR_filter` | `SOR > 10.0` | Strand odds ratio (relaxed for INDELs) |
| `QUAL_filter` | `QUAL < 30.0` | Raw call confidence |

> `ReadPosRankSum` and `MQ` filters are intentionally omitted for INDELs per GATK best practice.

---

## Key Output: multianno.csv Columns

| Column | Description |
|---|---|
| Chr, Start, End, Ref, Alt | Variant coordinates |
| Func.refGene | Gene function |
| Gene.refGene | Gene name |
| ExonicFunc.refGene | Exonic function (missense, nonsense, etc.) |
| AAChange.refGene | Amino acid change |
| gnomad41_exome | gnomAD exome allele frequency |
| clinvar_20250721 | ClinVar annotation |
| **Caller** | `BOTH`, `GATK_ONLY`, `DV_ONLY`, or `GATK` |
| **DeNovo** | `hiConfDeNovo:sample`, `loConfDeNovo:sample`, or `.` |
| **GT_sample** | `het`, `hom-alt`, `hom-ref`, or `missing` (one column per sample) |

---

## Notes

- **Lane merging:** When a sample has multiple lanes, all lane BAMs are merged with `samtools merge` before MarkDuplicates. Read groups distinguish lanes for BQSR. Supports 2+ lanes.
- **Multi-kit cohort safety:** `BCFTOOLS_SETGT` converts `DP=0` genotypes to `./. ` after joint calling. This prevents false de novo calls and allele frequency distortion when samples from different capture kits are jointly genotyped.
- **Dual caller concordance:** Variants called by both GATK and DeepVariant (`CALLER=BOTH`) have the highest confidence.
- **SpliceAI efficiency:** SpliceAI runs on the cohort-level VCF before per-sample splitting  scores are computed once per unique variant position.
- **Batch isolation:** Each batch's outputs are namespaced by `BATCH_NAME`. Running BATCH2 never overwrites BATCH1 results.

---

## Container Verification

Ensure these containers exist in the configured Singularity image directory:

```
fastqc_0.12.1.sif
fastp_0.23.4.sif
bwa-mem2-samtool2.3.1_1.22.1.sif
samtools_1.23.sif
gatk4_4.6.2.0.sif
deepvariant_1.9.0.sif
glnexus_1.2.7.sif
bcftool_1.23.sif
mosdepth_0.3.11.sif
multiqc_1.33.sif
peddy_0.4.8.sif
annovar_2020Jun08.sif
python_3.14.2.sif
spliceai.sif
```

### External Resources

For container build instructions, refer to the build repository:

[**BUILD_SINGULARITY_IMAGES**](https://github.com/Deep2106/BUILD_SINGULARITY_IMAGES/tree/main)

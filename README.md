
# WES/WGS Germline Variant Calling Pipeline v2.0

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
sample_id,fastq_1,fastq_2,family_id,paternal_id,maternal_id,sex,phenotype
child1,/path/to/child1_R1.fq.gz,/path/to/child1_R2.fq.gz,FAM001,father1,mother1,male,affected
father1,/path/to/father1_R1.fq.gz,/path/to/father1_R2.fq.gz,FAM001,0,0,male,unaffected
mother1,/path/to/mother1_R1.fq.gz,/path/to/mother1_R2.fq.gz,FAM001,0,0,female,unaffected
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
├── annotation/
│   ├── batch.hg38_multianno.csv    # ← Main output with GT, Caller, DeNovo columns
│   ├── batch.hg38_multianno.txt
│   └── batch_variants.vcf.gz        # Batch-filtered VCF
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

>>>>>
>>>>>

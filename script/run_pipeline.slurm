#!/bin/bash
#SBATCH --job-name=wes_v2_pipeline
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=168:00:00
#SBATCH --partition=ProdQ
#SBATCH --output=nextflow_%j.out
#SBATCH --error=nextflow_%j.err

# ============================================================================
# WES/WGS GERMLINE VARIANT CALLING PIPELINE v2.0 - RCSI HPC
# ============================================================================
# Features:
#   - Incremental GenomicsDB (auto-detection: CREATE vs UPDATE)
#   - Optional DeepVariant dual-caller with merged VCF
#   - Automatic PED file generation and management
#   - De novo mutation flagging (hiConfDeNovo/loConfDeNovo)
#   - Batch-specific annotation (only new samples)
#   - Automatic backup of GenomicsDB and PED before incremental runs
#
# Usage:
#   FIRST BATCH:  sbatch run_pipeline.slurm
#   NEXT BATCHES: Edit BATCH_NAME and INPUT_CSV, then sbatch run_pipeline.slurm
# ============================================================================

# ============================================================================
# USER CONFIGURATION - EDIT THESE FOR EACH BATCH
# ============================================================================

# Batch identifier (used for output organization)
BATCH_NAME="batch1"

# Input sample CSV for this batch
INPUT_CSV="sample_sheet_required.csv"

# Enable DeepVariant dual-caller? (true/false)
RUN_DEEPVARIANT="true"

# ============================================================================
# ENVIRONMENT SETUP
# ============================================================================

module purge
module load Java/17

# ============================================================================
# DIRECTORY SETUP
# ============================================================================

# Base project directory (on NFS)
PROJECT_DIR=/home/data/human_genetics/ONE4ALL

# Pipeline location
PIPELINE_DIR=${PROJECT_DIR}/TT/wes-pipeline2

# Nextflow work directory - MUST be on NFS (shared storage)
export NXF_WORK=${PROJECT_DIR}/work2
mkdir -p $NXF_WORK

# Output directory
OUTDIR=${PROJECT_DIR}/TEST_DVG/results_v2

# GenomicsDB path (persistent across batches)
GENOMICSDB_PATH=${PROJECT_DIR}/TEST_DVG/genomicsdb

# Master PED file path (auto-generated on first run, updated on subsequent runs)
MASTER_PED=${OUTDIR}/ped/master.ped

# ============================================================================
# SINGULARITY CONFIGURATION (per RCSI documentation)
# ============================================================================

# Library directory - pre-pulled SIF files
export NXF_SINGULARITY_LIBRARYDIR=/home/data/human_genetics/ONE4ALL/Docker_Images

# Cache directory for any additional downloads
export NXF_SINGULARITY_CACHEDIR=${PROJECT_DIR}/.singularity_cache
mkdir -p $NXF_SINGULARITY_CACHEDIR

# ============================================================================
# JAVA CONFIGURATION
# ============================================================================

# RCSI sets TMPDIR to local SSD, use it for Java temp files
export _JAVA_OPTIONS="-Djava.io.tmpdir=${TMPDIR:-/tmp}"

# ============================================================================
# NEXTFLOW OPTIONS
# ============================================================================

# Disable ANSI colors in logs (cleaner output)
export NXF_ANSI_LOG=false

# ============================================================================
# AUTO-DETECT MODE: INITIAL vs INCREMENTAL
# ============================================================================

# Check if GenomicsDB and master.ped exist
if [ -d "${GENOMICSDB_PATH}" ] && [ -f "${MASTER_PED}" ]; then
    RUN_MODE="INCREMENTAL"
    MASTER_PED_ARG="--master_ped ${MASTER_PED}"
    echo "MODE: INCREMENTAL (GenomicsDB and master.ped found)"
else
    RUN_MODE="INITIAL"
    MASTER_PED_ARG=""
    echo "MODE: INITIAL (First run - will create GenomicsDB and master.ped)"
fi

# ============================================================================
# BACKUP BEFORE INCREMENTAL RUN
# ============================================================================

if [ "${RUN_MODE}" == "INCREMENTAL" ]; then
    BACKUP_TIMESTAMP=$(date +%Y%m%d_%H%M%S)
    
    # Backup GenomicsDB
    if [ -d "${GENOMICSDB_PATH}" ]; then
        GENOMICSDB_BACKUP="${GENOMICSDB_PATH}_backup_${BACKUP_TIMESTAMP}"
        echo "Backing up GenomicsDB to: ${GENOMICSDB_BACKUP}"
        cp -rL "${GENOMICSDB_PATH}" "${GENOMICSDB_BACKUP}"
        if [ $? -eq 0 ]; then
            echo "GenomicsDB backup complete: $(du -sh ${GENOMICSDB_BACKUP} | cut -f1)"
        else
            echo "WARNING: GenomicsDB backup failed!"
        fi
    fi
    
    # Backup master.ped (additional backup in SLURM, module also creates one)
    if [ -f "${MASTER_PED}" ]; then
        PED_BACKUP="${MASTER_PED}_slurm_backup_${BACKUP_TIMESTAMP}"
        echo "Backing up master.ped to: ${PED_BACKUP}"
        cp "${MASTER_PED}" "${PED_BACKUP}"
    fi
    
    echo "------------------------------------------------------------------------------"
fi

# ============================================================================
# RUN PIPELINE
# ============================================================================

cd ${PIPELINE_DIR} || { echo "ERROR: Cannot cd to ${PIPELINE_DIR}"; exit 1; }

echo "=============================================================================="
echo "WES/WGS Germline Variant Calling Pipeline v2.0"
echo "=============================================================================="
echo "Date: $(date)"
echo "Node: $(hostname)"
echo "SLURM Job ID: ${SLURM_JOB_ID}"
echo "------------------------------------------------------------------------------"
echo "Batch Name: ${BATCH_NAME}"
echo "Input CSV: ${INPUT_CSV}"
echo "Run Mode: ${RUN_MODE}"
echo "DeepVariant: ${RUN_DEEPVARIANT}"
echo "------------------------------------------------------------------------------"
echo "TMPDIR (local SSD): ${TMPDIR}"
echo "NXF_WORK (NFS): ${NXF_WORK}"
echo "GENOMICSDB_PATH: ${GENOMICSDB_PATH}"
echo "MASTER_PED: ${MASTER_PED}"
echo "OUTPUT_DIR: ${OUTDIR}"
echo "------------------------------------------------------------------------------"
echo "NXF_SINGULARITY_LIBRARYDIR: ${NXF_SINGULARITY_LIBRARYDIR}"
echo "=============================================================================="

# Verify input CSV exists
if [ ! -f "${PROJECT_DIR}/${INPUT_CSV}" ]; then
    echo "ERROR: Input CSV not found: ${PROJECT_DIR}/${INPUT_CSV}"
    exit 1
fi

# List available SIF containers
echo "Available SIF containers:"
ls -la ${NXF_SINGULARITY_LIBRARYDIR}/*.sif 2>/dev/null | head -20 || echo "No .sif files found"
echo "=============================================================================="

# Run the pipeline
nextflow run main.nf \
    -profile slurm \
    -w ${NXF_WORK} \
    --input ${PROJECT_DIR}/${INPUT_CSV} \
    --fasta ${PROJECT_DIR}/reference/hg38.fasta \
    --fasta_fai ${PROJECT_DIR}/reference/hg38.fasta.fai \
    --fasta_dict ${PROJECT_DIR}/reference/hg38.dict \
    --bwa_index ${PROJECT_DIR}/reference \
    --targets_bed ${PROJECT_DIR}/reference/S000021_hg38_targets_extended.bed \
    --dbsnp ${PROJECT_DIR}/resources/hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf.gz \
    --dbsnp_tbi ${PROJECT_DIR}/resources/hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi \
    --known_indels ${PROJECT_DIR}/resources/hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz \
    --known_indels_tbi ${PROJECT_DIR}/resources/hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz.tbi \
    --annovar_db ${PROJECT_DIR}/Annovardb/humandb \
    --genomicsdb_path ${GENOMICSDB_PATH} \
    --outdir ${OUTDIR} \
    --run_deepvariant ${RUN_DEEPVARIANT} \
    ${MASTER_PED_ARG} \
    -resume

EXIT_CODE=$?

echo "=============================================================================="
echo "Pipeline finished with exit code: ${EXIT_CODE}"
echo "Date: $(date)"
echo "------------------------------------------------------------------------------"
if [ ${EXIT_CODE} -eq 0 ]; then
    echo "SUCCESS!"
    echo ""
    echo "Outputs:"
    echo "  - Annotated variants: ${OUTDIR}/annotation/"
    echo "  - Joint VCF: ${OUTDIR}/variants/joint/"
    echo "  - QC reports: ${OUTDIR}/qc/multiqc/"
    echo "  - Master PED: ${OUTDIR}/ped/master.ped"
    echo "  - GenomicsDB: ${GENOMICSDB_PATH}"
    echo ""
    echo "For next batch:"
    echo "  1. Edit BATCH_NAME and INPUT_CSV at top of this script"
    echo "  2. Run: sbatch run_pipeline.slurm"
else
    echo "FAILED - Check logs:"
    echo "  - Nextflow log: ${PIPELINE_DIR}/.nextflow.log"
    echo "  - SLURM output: nextflow_${SLURM_JOB_ID}.out"
    echo "  - SLURM error: nextflow_${SLURM_JOB_ID}.err"
fi
echo "=============================================================================="

exit ${EXIT_CODE}


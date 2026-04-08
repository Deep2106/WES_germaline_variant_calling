#!/bin/bash
#SBATCH --job-name=ciindi_b2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=168:00:00
#SBATCH --partition=ProdQ
#SBATCH --output=nextflow_%j.out
#SBATCH --error=nextflow_%j.err
#SBATCH --mail-user=deepakbharti@rcsi.com
#SBATCH --mail-type=END,FAIL


# ============================================================================
# WES/WGS GERMLINE VARIANT CALLING PIPELINE v2.0 - RCSI HPC
# ============================================================================
# Features:
#   - Incremental GenomicsDB (auto-detection: CREATE vs UPDATE)
#   - Optional DeepVariant dual-caller with merged VCF
#   - Automatic PED file generation and management
#   - De novo mutation flagging (hiConfDeNovo/loConfDeNovo)
#   - Batch-specific annotation (only new samples)
#   - Automatic backup of GenomicsDB, DeepVariant GVCFs, and PED before incremental runs
#   - Optional SpliceAI annotation
#
# Usage:
#   FIRST BATCH:  sbatch run_pipeline.slurm
#   NEXT BATCHES: Edit BATCH_NAME and INPUT_CSV, then sbatch run_pipeline.slurm
# ============================================================================

# This scrupts assumes /local/scratch is SSD hence every major operation done
# in /local/scratch. you may need to change configurations file if its not the
# case for your HPC

# ============================================================================
# USER CONFIGURATION - EDIT THESE FOR EACH BATCH
# ============================================================================

# Batch identifier (used for output organization and MultiQC reports)
BATCH_NAME="BATCH1"

# Input sample CSV for this batch
INPUT_CSV="SAMPLES.csv"

# Enable DeepVariant dual-caller? (true/false)
RUN_DEEPVARIANT="true"

# Enable SpliceAI annotation? (true/false)
RUN_SPLICEAI="true"

# ============================================================================
# ENVIRONMENT SETUP
# ============================================================================

module purge
module load Java/17

# ============================================================================
# DIRECTORY SETUP
# ============================================================================

# Base project directory (on NFS)
PROJECT_DIR=/path/to/project_directory    # 

# Pipeline location
PIPELINE_DIR=${PROJECT_DIR}/this/repository/location 

# Nextflow work directory - MUST be on NFS (shared storage)
export NXF_WORK=${PROJECT_DIR}/work ;
mkdir -p $NXF_WORK 

# Output directory
OUTDIR=${PROJECT_DIR}/results
mkdir -p $OUTDIR ;
# Derived paths - ALL consistent with OUTDIR
GENOMICSDB_PATH=${OUTDIR}/genomicsdb/genomicsdb 
DV_GVCF_PATH=${OUTDIR}/variants/deepvariant 
MASTER_PED=${OUTDIR}/ped/master.ped 

# Backup directory
BACKUP_BASE=${PROJECT_DIR}/Backups

# ============================================================================
# SINGULARITY CONFIGURATION (per RCSI documentation)
# ============================================================================

# Library directory - pre-pulled SIF files
export NXF_SINGULARITY_LIBRARYDIR=/path/to/docker/images  ## currently we are using prebuild singularity .sif containers #will add a short note how to make them #

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

# Check if GenomicsDB and master.ped exist if exists incrmenetal mode is automatically selected, in case first batch run is aborted or failed first delete the master_ped file
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
    BACKUP_DIR=${BACKUP_BASE}/${BATCH_NAME}_${BACKUP_TIMESTAMP}
    mkdir -p ${BACKUP_DIR}

    echo "=============================================================================="
    echo "CREATING BACKUPS BEFORE ${BATCH_NAME} RUN"
    echo "=============================================================================="

    # Backup GenomicsDB
    if [ -d "${GENOMICSDB_PATH}" ]; then
        echo "Backing up GenomicsDB..."
        cp -rL "${GENOMICSDB_PATH}" "${BACKUP_DIR}/genomicsdb"
        if [ $? -eq 0 ]; then
            echo "  GenomicsDB backup: ${BACKUP_DIR}/genomicsdb ($(du -sh ${BACKUP_DIR}/genomicsdb | cut -f1))"
        else
            echo "  WARNING: GenomicsDB backup failed!"
        fi
    fi

    # Backup DeepVariant GVCFs
    if [ -d "${DV_GVCF_PATH}" ]; then
        echo "Backing up DeepVariant GVCFs..."
        cp -rL "${DV_GVCF_PATH}" "${BACKUP_DIR}/deepvariant_gvcfs"
        if [ $? -eq 0 ]; then
            DV_COUNT=$(ls ${BACKUP_DIR}/deepvariant_gvcfs/*.g.vcf.gz 2>/dev/null | wc -l)
            echo "  DeepVariant GVCFs backup: ${BACKUP_DIR}/deepvariant_gvcfs (${DV_COUNT} files)"
        else
            echo "  WARNING: DeepVariant GVCFs backup failed!"
        fi
    fi

    # Backup master.ped
    if [ -f "${MASTER_PED}" ]; then
        cp "${MASTER_PED}" "${BACKUP_DIR}/master.ped"
        echo "  master.ped backup: ${BACKUP_DIR}/master.ped"
    fi

    # Log sample counts
    if [ -f "${GENOMICSDB_PATH}/callset.json" ]; then
        SAMPLE_COUNT=$(grep "row_idx" "${GENOMICSDB_PATH}/callset.json" | wc -l)
        echo "  Current samples in GenomicsDB: ${SAMPLE_COUNT}"
    fi

    echo "=============================================================================="
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
echo "SpliceAI: ${RUN_SPLICEAI}"
echo "------------------------------------------------------------------------------"
echo "TMPDIR (local SSD): ${TMPDIR}"
echo "NXF_WORK (NFS): ${NXF_WORK}"
echo "OUTPUT_DIR: ${OUTDIR}"
echo "GENOMICSDB_PATH: ${GENOMICSDB_PATH}"
echo "DV_GVCF_PATH: ${DV_GVCF_PATH}"
echo "MASTER_PED: ${MASTER_PED}"
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

# its better to create paded bed file for the run eg 30 BP up/down 
# you may choose your own truth set 
# Run the pipeline
nextflow run main.nf \
    -profile slurm \
    -w ${NXF_WORK} \
    --input ${PROJECT_DIR}/${INPUT_CSV} \
    --batch_name ${BATCH_NAME} \
    --fasta /path/to/reference/hg38.fasta \
    --fasta_fai /path/to/reference/hg38.fasta.fai \
    --fasta_dict /path/to/reference/hg38.dict \
    --bwa_index /path/to/reference \
    --targets_bed /path/to/bed_file \  
    --dbsnp /path/to/dbsnpxxx.vcf.gz \
    --dbsnp_tbi /path/to/dbsnpxxx.vcf.gz.tbi \
    --known_indels /path/to/resources/hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz \
    --known_indels_tbi /path/to/resources/hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz.tbi \
    --annovar_db /path/to/Annovardb/humandb \
    --genomicsdb_path ${GENOMICSDB_PATH} \
    --dv_gvcf_path ${DV_GVCF_PATH} \
    --outdir ${OUTDIR} \
    --run_deepvariant ${RUN_DEEPVARIANT} \
    --run_spliceai ${RUN_SPLICEAI} \
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
    echo "  - QC reports: ${OUTDIR}/qc/multiqc/${BATCH_NAME}/"
    echo "  - Master PED: ${MASTER_PED}"
    echo "  - GenomicsDB: ${GENOMICSDB_PATH}"
    echo "  - DeepVariant GVCFs: ${DV_GVCF_PATH}"
    if [ -n "${BACKUP_DIR}" ]; then
        echo "  - Pre-run backup: ${BACKUP_DIR}"
    fi
    echo ""

    # Show updated sample count
    if [ -f "${GENOMICSDB_PATH}/callset.json" ]; then
        NEW_SAMPLE_COUNT=$(grep "row_idx" "${GENOMICSDB_PATH}/callset.json" | wc -l)
        echo "Total samples in GenomicsDB after ${BATCH_NAME}: ${NEW_SAMPLE_COUNT}"
    fi
    
    # Show DeepVariant GVCF count
    if [ -d "${DV_GVCF_PATH}" ]; then
        DV_GVCF_COUNT=$(ls ${DV_GVCF_PATH}/*.g.vcf.gz 2>/dev/null | wc -l)
        echo "Total DeepVariant GVCFs: ${DV_GVCF_COUNT}"
    fi
    
    echo ""
    echo "For next batch:"
    echo "  1. Edit BATCH_NAME and INPUT_CSV at top of this script"
    echo "  2. Run: sbatch run_pipeline.slurm"
else
    echo "FAILED - Check logs:"
    echo "  - Nextflow log: ${PIPELINE_DIR}/.nextflow.log"
    echo "  - SLURM output: nextflow_${SLURM_JOB_ID}.out"
    echo "  - SLURM error: nextflow_${SLURM_JOB_ID}.err"
    if [ -n "${BACKUP_DIR}" ]; then
        echo ""
        echo "To restore from backup:"
        echo "  rm -rf ${GENOMICSDB_PATH}"
        echo "  cp -r ${BACKUP_DIR}/genomicsdb ${GENOMICSDB_PATH}"
        echo "  cp -r ${BACKUP_DIR}/deepvariant_gvcfs ${DV_GVCF_PATH}"
        echo "  cp ${BACKUP_DIR}/master.ped ${MASTER_PED}"
    fi
fi
echo "=============================================================================="

exit ${EXIT_CODE}

## clean your work directory

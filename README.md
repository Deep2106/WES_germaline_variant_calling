# WES_germaline_variant_calling
Containerized Nextflow pipeline for WES, easily extendable to WGS (via VQSR). It integrates FastQC/MultiQC for QC, fastp for preprocessing, BWA-MEM2 for alignment, GATK 4.2.6 for variant calling, and ANNOVAR for annotation, with optional DeepVariant and SpliceAI support.

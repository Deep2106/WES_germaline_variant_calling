/*
========================================================================================
    GATK GENOMICSDBIMPORT MODULE
========================================================================================
    FIXES:
    - Remove symlink before copying back (publishDir needs real directory)
    - Replaced broken sed callset.json parser with Python JSON parser
      Reason: GATK >= 4.5 / GenomicsDB >= 1.4 writes array format:
        {"callsets": [{"sample_name": "...", "row_idx": "0", ...}]}
      Old sed assumed dict format: {"SAMPLE": {"row_idx": 0}}
      Mismatch caused 0 existing samples detected -> all flagged as NEW ->
      GenomicsDBException: Duplicate sample on UPDATE runs
    - Parser handles both array and dict formats with hard fail on parse
      error (no silent suppression via || true)
========================================================================================
*/

process GATK_GENOMICSDBIMPORT {
    tag "genomicsdb_${update_db ? 'UPDATE' : 'CREATE'}"
    label 'process_high'
    container "${params.containers.gatk4}"
    
    publishDir "${params.outdir}/genomicsdb", mode: 'copy', overwrite: true
    
    input:
    path gvcfs
    path gvcf_indices
    path intervals
    path existing_db
    val update_db
    
    output:
    path "genomicsdb", emit: genomicsdb
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def avail_mem = task.memory ? task.memory.toGiga() - 8 : 56
    def java_opts = "-Xmx${avail_mem}g -XX:+UseParallelGC -XX:ParallelGCThreads=4 -Djava.io.tmpdir=\$TMPDIR"
    def intervals_arg = intervals.name != 'NO_FILE' ? "--intervals ${intervals}" : ""
    def reader_threads = Math.max(1, task.cpus - 2)
    def parallel_intervals = Math.max(1, (task.cpus / 4).intValue())
    
    if (update_db && existing_db.name != 'NO_DB') {
        """
        #!/bin/bash
        set -euo pipefail
        
        if [ -n "\${SLURM_JOB_ID:-}" ] && [ -d "/local/scratch/\${SLURM_JOB_ID}" ]; then
            export TMPDIR="/local/scratch/\${SLURM_JOB_ID}"
            LOCAL_SCRATCH="/local/scratch/\${SLURM_JOB_ID}/gdb_update"
        elif [ -n "\${TMPDIR:-}" ] && [ -d "\${TMPDIR}" ]; then
            LOCAL_SCRATCH="\${TMPDIR}/gdb_update"
        else
            export TMPDIR=\$(mktemp -d)
            LOCAL_SCRATCH="\${TMPDIR}/gdb_update"
        fi
        mkdir -p "\${LOCAL_SCRATCH}"
        
        echo "=============================================="
        echo "GENOMICSDB UPDATE MODE"
        echo "=============================================="
        echo "Existing DB: ${existing_db}"
        
        # Check if input is symlink (Nextflow staging)
        if [ -L "${existing_db}" ]; then
            echo "Input is symlink -> \$(readlink -f "${existing_db}")"
        fi
        
        echo ""
        echo "Checking for existing samples in GenomicsDB..."

        # Extract sample names from callset.json
        # Supports both formats written by different GenomicsDB versions:
        #   Array format  (GenomicsDB >= 1.4 / GATK >= 4.5):
        #     {"callsets": [{"sample_name": "SAMPLE", "row_idx": "0", ...}, ...]}
        #   Dict format   (GenomicsDB < 1.4 / GATK < 4.5):
        #     {"callsets": {"SAMPLE": {"row_idx": 0}, ...}}
        if [ -f "${existing_db}/callset.json" ]; then
            python3 - "${existing_db}/callset.json" > existing_samples.txt <<'PYEOF'
import json, sys

callset_file = sys.argv[1]
try:
    with open(callset_file) as fh:
        data = json.load(fh)
except Exception as e:
    print(f"ERROR: Failed to parse {callset_file}: {e}", file=sys.stderr)
    sys.exit(1)

callsets = data.get("callsets", [])

if isinstance(callsets, list):
    # Array format: [{"sample_name": "...", ...}, ...]
    for entry in callsets:
        name = entry.get("sample_name", "").strip()
        if name:
            print(name)
elif isinstance(callsets, dict):
    # Dict format: {"SAMPLE": {"row_idx": ...}, ...}
    for name in callsets:
        name = name.strip()
        if name:
            print(name)
else:
    print(f"ERROR: Unexpected callsets type: {type(callsets)}", file=sys.stderr)
    sys.exit(1)
PYEOF

            if [ \$? -ne 0 ]; then
                echo "ERROR: callset.json parsing failed - cannot safely determine existing samples" >&2
                echo "       Refusing to continue to prevent duplicate-sample corruption of GenomicsDB" >&2
                exit 1
            fi

            EXISTING_COUNT=\$(wc -l < existing_samples.txt)
            echo "Found \${EXISTING_COUNT} existing samples (JSON parsed):"
            cat existing_samples.txt
        else
            touch existing_samples.txt
            EXISTING_COUNT=0
            echo "WARNING: No callset.json found in ${existing_db} - treating DB as empty"
        fi
        
        echo ""
        echo "Building sample map..."
        > sample_map.txt
        SKIPPED=0
        ADDED=0
        
        for gvcf in *.g.vcf.gz; do
            if [ -f "\$gvcf" ]; then
                sample_name=\$(bcftools query -l "\$gvcf" | head -1)
                if grep -qx "\${sample_name}" existing_samples.txt; then
                    echo "  SKIP: \${sample_name} (already in DB)"
                    ((SKIPPED++)) || true
                else
                    printf "%s\\t%s\\n" "\${sample_name}" "\$(readlink -f "\$gvcf")" >> sample_map.txt
                    echo "  ADD: \${sample_name}"
                    ((ADDED++)) || true
                fi
            fi
        done
        
        echo ""
        echo "Summary: New=\${ADDED}, Skipped=\${SKIPPED}"
        
        if [ \${ADDED} -eq 0 ]; then
            echo "No new samples to add"
            # Remove symlink if exists, then copy
            [ -L genomicsdb ] && rm genomicsdb
            [ -d genomicsdb ] && rm -rf genomicsdb
            cp -rL "${existing_db}" genomicsdb
            echo '"${task.process}":' > versions.yml
            echo '    gatk4: "4.6.2.0"' >> versions.yml
            echo '    mode: "update_skipped"' >> versions.yml
            exit 0
        fi
        
        LOCAL_GENOMICSDB="\${LOCAL_SCRATCH}/genomicsdb_temp"
        echo ""
        echo "Copying existing DB to local scratch..."
        cp -rL "${existing_db}" "\${LOCAL_GENOMICSDB}"
        echo "Local copy size: \$(du -sh "\${LOCAL_GENOMICSDB}" | cut -f1)"
        
        echo ""
        echo "Running GenomicsDBImport UPDATE..."
        gatk --java-options "${java_opts}" GenomicsDBImport \\
            --sample-name-map sample_map.txt \\
            --genomicsdb-update-workspace-path "\${LOCAL_GENOMICSDB}" \\
            ${intervals_arg} \\
            --tmp-dir "\$TMPDIR" \\
            --merge-input-intervals true \\
            --reader-threads ${reader_threads} \\
            --max-num-intervals-to-import-in-parallel ${parallel_intervals} \\
            --genomicsdb-shared-posixfs-optimizations true \\
            --bypass-feature-reader true \\
            --batch-size 50
        
        echo ""
        echo "Copying updated DB back to work directory..."
        # CRITICAL: Remove symlink first, otherwise cp goes to symlink target
        [ -L genomicsdb ] && rm genomicsdb
        [ -d genomicsdb ] && rm -rf genomicsdb
        cp -rL "\${LOCAL_GENOMICSDB}" genomicsdb
        rm -rf "\${LOCAL_SCRATCH}"
        
        echo ""
        echo "UPDATE complete: \${ADDED} samples added"
        echo "Final GenomicsDB size: \$(du -sh genomicsdb | cut -f1)"
        
        # Verify it's a real directory, not symlink
        if [ -L genomicsdb ]; then
            echo "ERROR: genomicsdb is still a symlink!" >&2
            exit 1
        fi
        
        echo '"${task.process}":' > versions.yml
        echo "    gatk4: \\"\$(gatk --version 2>&1 | head -n1 | sed 's/^.*(GATK) v//' | sed 's/ .*\$//')\\"" >> versions.yml
        echo '    mode: "update"' >> versions.yml
        echo "    samples_added: \${ADDED}" >> versions.yml
        """
    } else {
        """
        #!/bin/bash
        set -euo pipefail
        
        if [ -n "\${SLURM_JOB_ID:-}" ] && [ -d "/local/scratch/\${SLURM_JOB_ID}" ]; then
            export TMPDIR="/local/scratch/\${SLURM_JOB_ID}"
            LOCAL_SCRATCH="/local/scratch/\${SLURM_JOB_ID}/gdb_create"
        elif [ -n "\${TMPDIR:-}" ] && [ -d "\${TMPDIR}" ]; then
            LOCAL_SCRATCH="\${TMPDIR}/gdb_create"
        else
            export TMPDIR=\$(mktemp -d)
            LOCAL_SCRATCH="\${TMPDIR}/gdb_create"
        fi
        mkdir -p "\${LOCAL_SCRATCH}"
        
        echo "=============================================="
        echo "GENOMICSDB CREATE MODE"
        echo "=============================================="
        
        # Clean up any existing genomicsdb (symlink or directory)
        [ -L genomicsdb ] && rm genomicsdb
        [ -d genomicsdb ] && rm -rf genomicsdb
        
        > sample_map.txt
        for gvcf in *.g.vcf.gz; do
            if [ -f "\$gvcf" ]; then
                sample_name=\$(bcftools query -l "\$gvcf" | head -1)
                printf "%s\\t%s\\n" "\${sample_name}" "\$(readlink -f "\$gvcf")" >> sample_map.txt
            fi
        done
        
        echo "Samples:"
        cat sample_map.txt
        SAMPLE_COUNT=\$(wc -l < sample_map.txt)
        
        LOCAL_GENOMICSDB="\${LOCAL_SCRATCH}/genomicsdb_temp"
        
        echo ""
        echo "Running GenomicsDBImport CREATE..."
        gatk --java-options "${java_opts}" GenomicsDBImport \\
            --sample-name-map sample_map.txt \\
            --genomicsdb-workspace-path "\${LOCAL_GENOMICSDB}" \\
            ${intervals_arg} \\
            --tmp-dir "\$TMPDIR" \\
            --merge-input-intervals true \\
            --reader-threads ${reader_threads} \\
            --max-num-intervals-to-import-in-parallel ${parallel_intervals} \\
            --genomicsdb-shared-posixfs-optimizations true \\
            --bypass-feature-reader true \\
            --batch-size 50
        
        echo ""
        echo "Copying DB to work directory..."
        cp -rL "\${LOCAL_GENOMICSDB}" genomicsdb
        rm -rf "\${LOCAL_SCRATCH}"
        
        echo ""
        echo "CREATE complete: \${SAMPLE_COUNT} samples"
        echo "Final GenomicsDB size: \$(du -sh genomicsdb | cut -f1)"
        
        echo '"${task.process}":' > versions.yml
        echo "    gatk4: \\"\$(gatk --version 2>&1 | head -n1 | sed 's/^.*(GATK) v//' | sed 's/ .*\$//')\\"" >> versions.yml
        echo '    mode: "create"' >> versions.yml
        echo "    samples_imported: \${SAMPLE_COUNT}" >> versions.yml
        """
    }
}

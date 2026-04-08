/*
========================================================================================
    ANNOVAR MODULE - WITH GENOTYPE (GT) + SPLICEAI INFORMATION
========================================================================================
    Functional annotation of genetic variants using ANNOVAR

    Annotates variants with:
    - Gene-based annotation (refGene)
    - Region-based annotation (cytoBand)
    - Filter-based annotation (gnomAD, ClinVar, dbNSFP, etc.)
    - SpliceAI scores extracted from VCF INFO field

    Input:
        meta        - Sample metadata [sample_id, family_id, etc.]
        vcf         - VCF file to annotate (single sample, merged SNP+INDEL)
        tbi         - VCF index
        annovar_db  - ANNOVAR database directory

    Output:
        txt         - Annotated variants (multianno.txt)
        vcf         - Annotated VCF
        csv         - Annotated variants CSV (with GT, SpliceAI columns)
        avinput     - ANNOVAR input file with GT information
        versions    - Software versions

    Note:
        - .avinput contains GT (het/hom) in column 6 when using -withzyg
        - .csv output includes GT, DeNovo, Caller, and SpliceAI columns
        - SpliceAI scores: DS_AG, DS_AL, DS_DG, DS_DL, max score, gene
        - To filter in R:
            SNPs:   filter(nchar(Ref)==1 & nchar(Alt)==1 & Ref!="-" & Alt!="-")
            INDELs: filter(nchar(Ref)!=1 | nchar(Alt)!=1 | Ref=="-" | Alt=="-")
----------------------------------------------------------------------------------------
*/

process ANNOVAR {
    tag "${meta.sample_id}"
    label 'process_medium'

    container "${params.containers?.annovar ?: 'biocontainers/annovar:2020Jun08'}"

    publishDir "${params.outdir}/annotation/per_sample",
        mode: params.publish_dir_mode ?: 'copy',
        pattern: "*.{csv,avinput}", failOnError: false

    input:
    tuple val(meta), path(vcf), path(tbi)
    path annovar_db

    output:
    tuple val(meta), path("*.hg38_multianno.txt") , emit: txt
    tuple val(meta), path("*.hg38_multianno.vcf") , emit: vcf, optional: true
    tuple val(meta), path("*.hg38_multianno.csv") , emit: csv
    tuple val(meta), path("*.avinput")            , emit: avinput
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = "${meta.sample_id}"

    // ANNOVAR parameters from config - with validated defaults
    def buildver = params.annovar_buildver ?: 'hg38'
    def protocol = params.annovar_protocol ?: 'refGene,cytoBand,genomicSuperDups,exac03,gnomad41_exome,gnomad41_genome,avsnp151,dbnsfp47a,dbscsnv11,mcap,clinvar_20250721,revel,gene4denovo201907'
    def operation = params.annovar_operation ?: 'g,r,r,f,f,f,f,f,f,f,f,f,f'
    def argument = params.annovar_argument ?: "'-splicing 5',,,,,,,,,,,,"

    // Count protocols and validate
    def proto_count = protocol.split(',').size()
    def op_count = operation.split(',').size()

    """
    #!/bin/bash
    set -euo pipefail
    export PATH="/home/TOOLS/tools/annovar/current/bin:\$PATH"
    echo "=== ANNOVAR Annotation with GT + SpliceAI Information ==="
    echo "Sample: ${prefix}"
    echo "Input VCF: ${vcf}"
    echo "Build: ${buildver}"
    echo "Protocols (${proto_count}): ${protocol}"
    echo "Operations (${op_count}): ${operation}"

    # Validate protocol/operation count match
    if [ ${proto_count} -ne ${op_count} ]; then
        echo "ERROR: Protocol count (${proto_count}) != Operation count (${op_count})"
        exit 1
    fi

    # =========================================================================
    # Step 1: Extract SpliceAI scores from VCF INFO field
    # SpliceAI format: SpliceAI=A|GENE|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL
    # Creates a lookup: chr_pos_ref_alt -> scores
    # =========================================================================
    echo ""
    echo "=== Step 1: Extracting SpliceAI scores from VCF ==="

    zcat ${vcf} | grep -v "^#" | awk -F'\\t' '{
        chr = \$1; pos = \$2; ref = \$4; alt = \$5
        spliceai = "."
        n = split(\$8, info_fields, ";")
        for (i = 1; i <= n; i++) {
            if (index(info_fields[i], "SpliceAI=") == 1) {
                spliceai = substr(info_fields[i], 10)
                break
            }
        }
        if (spliceai != ".") {
            n2 = split(spliceai, sp, "|")
            if (n2 >= 6) {
                gene   = sp[2]
                ds_ag  = sp[3]; ds_al = sp[4]; ds_dg = sp[5]; ds_dl = sp[6]
                max_ds = ds_ag + 0
                if (ds_al + 0 > max_ds) max_ds = ds_al + 0
                if (ds_dg + 0 > max_ds) max_ds = ds_dg + 0
                if (ds_dl + 0 > max_ds) max_ds = ds_dl + 0
                key = chr "_" pos "_" pos "_" ref "_" alt
                print key "\\t" gene "\\t" ds_ag "\\t" ds_al "\\t" ds_dg "\\t" ds_dl "\\t" max_ds "\\t" spliceai
            }
        }
    }' > spliceai_lookup.txt

    SPLICEAI_COUNT=\$(wc -l < spliceai_lookup.txt)
    echo "Extracted SpliceAI scores for \${SPLICEAI_COUNT} variants"

    if [ "\${SPLICEAI_COUNT}" -gt 0 ]; then
        echo "Sample SpliceAI entries:"
        head -3 spliceai_lookup.txt
    fi

    # =========================================================================
    # Step 2: Convert VCF to ANNOVAR input format WITH ZYGOSITY (-withzyg)
    # =========================================================================
    echo ""
    echo "=== Step 2: Converting VCF to ANNOVAR format with zygosity ==="

    convert2annovar.pl \\
        -format vcf4 \\
        -includeinfo \\
        -withzyg \\
        ${vcf} \\
        > ${prefix}.avinput

    if [ ! -s ${prefix}.avinput ]; then
        echo "WARNING: Empty or missing .avinput file"
        touch ${prefix}.avinput
    else
        echo "Created ${prefix}.avinput with \$(wc -l < ${prefix}.avinput) variants"
        echo "Sample .avinput content (first 3 lines):"
        head -3 ${prefix}.avinput || true
    fi

    # Save a copy of avinput
    cp ${prefix}.avinput ${prefix}.avinput.original

    # =========================================================================
    # Step 3: Run ANNOVAR table annotation
    # =========================================================================
    echo ""
    echo "=== Step 3: Running ANNOVAR table annotation ==="

    table_annovar.pl \\
        ${prefix}.avinput \\
        ${annovar_db} \\
        -buildver ${buildver} \\
        -out ${prefix} \\
        -protocol ${protocol} \\
        -operation ${operation} \\
        -argument ${argument} \\
        -nastring . \\
        -polish \\
        -otherinfo \\
        ${args}

    # Restore original avinput if modified
    if [ ! -s ${prefix}.avinput ] && [ -s ${prefix}.avinput.original ]; then
        echo "Restoring .avinput from backup"
        mv ${prefix}.avinput.original ${prefix}.avinput
    else
        rm -f ${prefix}.avinput.original
    fi

    # =========================================================================
    # Step 4: Generate CSV with GT + SpliceAI columns using Perl
    # Perl is guaranteed available (ANNOVAR requires it)
    # Avoids all Nextflow/bash dollar-sign escaping issues
    # =========================================================================
    echo ""
    echo "=== Step 4: Generating CSV with GT and SpliceAI information ==="

    if [ ! -f ${prefix}.${buildver}_multianno.txt ]; then
        echo "ERROR: ${prefix}.${buildver}_multianno.txt not found!"
        ls -la ${prefix}.* || true
        exit 1
    fi

    # Create GT lookup from avinput
    awk -F'\\t' '{
        key = \$1"_"\$2"_"\$3"_"\$4"_"\$5
        gt = \$6
        print key"\\t"gt
    }' ${prefix}.avinput > gt_lookup.txt

    echo "Created GT lookup with \$(wc -l < gt_lookup.txt) entries"

    # Run Perl script for CSV generation
    perl -e '
    use strict;
    use warnings;

    # Load GT lookup
    my %gt_map;
    open(my \$gt_fh, "<", "gt_lookup.txt") or die "Cannot open gt_lookup.txt";
    while (<\$gt_fh>) {
        chomp;
        my (\$key, \$gt) = split(/\\t/, \$_, 2);
        \$gt_map{\$key} = \$gt if defined \$key && defined \$gt;
    }
    close(\$gt_fh);

    # Load SpliceAI lookup
    my (%sp_gene, %sp_ag, %sp_al, %sp_dg, %sp_dl, %sp_max, %sp_raw);
    if (-f "spliceai_lookup.txt" && -s "spliceai_lookup.txt") {
        open(my \$sp_fh, "<", "spliceai_lookup.txt") or die "Cannot open spliceai_lookup.txt";
        while (<\$sp_fh>) {
            chomp;
            my @f = split(/\\t/, \$_);
            next unless scalar @f >= 8;
            \$sp_gene{\$f[0]} = \$f[1];
            \$sp_ag{\$f[0]}   = \$f[2];
            \$sp_al{\$f[0]}   = \$f[3];
            \$sp_dg{\$f[0]}   = \$f[4];
            \$sp_dl{\$f[0]}   = \$f[5];
            \$sp_max{\$f[0]}  = \$f[6];
            \$sp_raw{\$f[0]}  = \$f[7];
        }
        close(\$sp_fh);
    }

    # Counters
    my (\$gt_matched, \$gt_unmatched, \$sp_matched, \$sp_unmatched) = (0, 0, 0, 0);

    # Find Otherinfo column index
    my \$otherinfo_col = -1;

    # Process multianno.txt
    my \$input_file = shift @ARGV;
    my \$output_file = shift @ARGV;
    open(my \$in, "<", \$input_file) or die "Cannot open \$input_file";
    open(my \$out, ">", \$output_file) or die "Cannot open \$output_file";

    my \$line_num = 0;
    while (<\$in>) {
        chomp;
        \$line_num++;
        my @cols = split(/\\t/, \$_, -1);

        if (\$line_num == 1) {
            # Find Otherinfo column
            for my \$i (0 .. \$#cols) {
                if (\$cols[\$i] =~ /Otherinfo/) {
                    \$otherinfo_col = \$i;
                    last;
                }
            }
            my \$last_col = (\$otherinfo_col > 0) ? \$otherinfo_col - 1 : \$#cols;

            # Print header
            my @hdr;
            for my \$i (0 .. \$last_col) {
                push @hdr, csv_quote(\$cols[\$i]);
            }
            # Insert GT, DeNovo, Caller after column 5 (index 4)
            splice(@hdr, 5, 0, "GT", "DeNovo", "Caller");
            # Append SpliceAI columns
            push @hdr, "SpliceAI_gene", "SpliceAI_DS_AG", "SpliceAI_DS_AL",
                       "SpliceAI_DS_DG", "SpliceAI_DS_DL", "SpliceAI_max", "SpliceAI_raw";
            print \$out join(",", @hdr) . "\\n";
            next;
        }

        # Data rows
        my \$last_col = (\$otherinfo_col > 0) ? \$otherinfo_col - 1 : \$#cols;
        my \$key = join("_", @cols[0..4]);

        # GT lookup
        my \$gt = exists \$gt_map{\$key} ? \$gt_map{\$key} : ".";
        if (\$gt ne ".") { \$gt_matched++ } else { \$gt_unmatched++ }

        # De novo and caller from Otherinfo columns
        my \$dn = ".";
        my \$cl = "GATK";
        my \$full_line = \$_;
        if (\$full_line =~ /hiConfDeNovo/)  { \$dn = "hiConfDeNovo" }
        elsif (\$full_line =~ /loConfDeNovo/) { \$dn = "loConfDeNovo" }
        if (\$full_line =~ /CALLER=BOTH/)      { \$cl = "BOTH" }
        elsif (\$full_line =~ /CALLER=GATK_ONLY/) { \$cl = "GATK_ONLY" }
        elsif (\$full_line =~ /CALLER=DV_ONLY/)   { \$cl = "DV_ONLY" }

        # SpliceAI lookup
        my \$sg = exists \$sp_gene{\$key} ? \$sp_gene{\$key} : ".";
        my \$sa = exists \$sp_ag{\$key}   ? \$sp_ag{\$key}   : ".";
        my \$sl = exists \$sp_al{\$key}   ? \$sp_al{\$key}   : ".";
        my \$sd = exists \$sp_dg{\$key}   ? \$sp_dg{\$key}   : ".";
        my \$se = exists \$sp_dl{\$key}   ? \$sp_dl{\$key}   : ".";
        my \$sm = exists \$sp_max{\$key}  ? \$sp_max{\$key}  : ".";
        my \$sr = exists \$sp_raw{\$key}  ? \$sp_raw{\$key}  : ".";
        if (\$sg ne ".") { \$sp_matched++ } else { \$sp_unmatched++ }

        # Build output row
        my @row;
        for my \$i (0 .. \$last_col) {
            push @row, csv_quote(\$cols[\$i] // ".");
        }
        # Insert GT, DeNovo, Caller after column 5
        splice(@row, 5, 0, csv_quote(\$gt), csv_quote(\$dn), csv_quote(\$cl));
        # Append SpliceAI
        push @row, csv_quote(\$sg), csv_quote(\$sa), csv_quote(\$sl),
                   csv_quote(\$sd), csv_quote(\$se), csv_quote(\$sm), csv_quote(\$sr);

        print \$out join(",", @row) . "\\n";
    }
    close(\$in);
    close(\$out);

    # Print stats to stderr
    print STDERR "GT matching: \$gt_matched matched, \$gt_unmatched unmatched\\n";
    print STDERR "SpliceAI matching: \$sp_matched matched, \$sp_unmatched unmatched\\n";

    sub csv_quote {
        my \$f = shift // ".";
        if (\$f =~ /,/) {
            \$f =~ s/"/""/g;
            return "\\"" . \$f . "\\"";
        }
        return \$f;
    }
    ' ${prefix}.${buildver}_multianno.txt ${prefix}.${buildver}_multianno.csv

    # Cleanup
    rm -f gt_lookup.txt spliceai_lookup.txt

    # =========================================================================
    # Summary statistics
    # =========================================================================
    echo ""
    echo "=== Summary ==="

    # Count de novo variants (column 7 = DeNovo after header)
    denovo_count=\$(awk -F',' 'NR>1 && \$7!="." {count++} END {print count+0}' ${prefix}.${buildver}_multianno.csv)
    echo "De novo variants found: \${denovo_count}"

    # Count caller distribution (column 8 = Caller)
    echo "Caller distribution:"
    awk -F',' 'NR>1 {caller[\$8]++} END {for (c in caller) print "  " c ": " caller[c]}' ${prefix}.${buildver}_multianno.csv

    # Count SpliceAI annotated variants by checking SpliceAI_max (second to last column)
    spliceai_annotated=\$(awk -F',' 'NR>1 { n=split(\$0,a,","); if (a[n-1] != ".") count++ } END {print count+0}' ${prefix}.${buildver}_multianno.csv)
    echo "Variants with SpliceAI scores: \${spliceai_annotated}"

    # Count high-impact splice variants
    high_impact_splice=\$(awk -F',' 'NR>1 { n=split(\$0,a,","); if (a[n-1]+0 >= 0.2) count++ } END {print count+0}' ${prefix}.${buildver}_multianno.csv)
    echo "High-impact splice variants (max_DS >= 0.2): \${high_impact_splice}"

    # Column counts
    txt_cols=\$(head -1 ${prefix}.${buildver}_multianno.txt | tr '\\t' '\\n' | wc -l)
    csv_cols=\$(head -1 ${prefix}.${buildver}_multianno.csv | tr ',' '\\n' | wc -l)
    echo "Columns: .txt has \${txt_cols}, .csv has \${csv_cols} (GT/DeNovo/Caller/SpliceAI added, Otherinfo removed)"

    # Show CSV header
    echo ""
    echo "CSV columns:"
    head -1 ${prefix}.${buildver}_multianno.csv | tr ',' '\\n' | cat -n || true

    # Final output
    echo ""
    echo "=== ANNOVAR + SpliceAI Complete ==="
    echo "Output files:"
    ls -lh ${prefix}.${buildver}_multianno.* ${prefix}.avinput 2>/dev/null || true

    total_txt=\$(tail -n +2 ${prefix}.${buildver}_multianno.txt | wc -l || echo "0")
    total_csv=\$(tail -n +2 ${prefix}.${buildver}_multianno.csv | wc -l || echo "0")
    total_avinput=\$(wc -l < ${prefix}.avinput || echo "0")

    echo ""
    echo "Variant counts:"
    echo "  .txt: \${total_txt} variants"
    echo "  .csv: \${total_csv} variants"
    echo "  .avinput: \${total_avinput} variants"

    # Versions
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        annovar: \$(perl -e 'print "2020Jun08"')
        spliceai_integrated: true
    END_VERSIONS

    # Explicit successful exit
    exit 0
    """

    stub:
    def prefix = "${meta.sample_id}"
    def buildver = params.annovar_buildver ?: 'hg38'
    """
    touch ${prefix}.${buildver}_multianno.txt
    touch ${prefix}.${buildver}_multianno.vcf
    echo "Chr,Start,End,Ref,Alt,GT,DeNovo,Caller,Gene,SpliceAI_gene,SpliceAI_DS_AG,SpliceAI_DS_AL,SpliceAI_DS_DG,SpliceAI_DS_DL,SpliceAI_max,SpliceAI_raw" > ${prefix}.${buildver}_multianno.csv
    echo -e "chr1\\t100\\t100\\tA\\tG\\thet" > ${prefix}.avinput

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        annovar: stub
        spliceai_integrated: true
    END_VERSIONS
    """
}

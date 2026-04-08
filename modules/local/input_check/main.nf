/*
    INPUT CHECK MODULE
    Validates sample CSV and generates/merges PED files
    
    CSV format (standard PED order):
    family_id,sample_id,paternal_id,maternal_id,sex,phenotype,fastq_r1,fastq_r2
    
    BACKUP: Creates timestamped backup of master.ped before merging
*/

process INPUT_CHECK {
    tag "input_check"
    label 'process_low'
    container "${params.containers.python}"
    
    publishDir "${params.outdir}/ped", mode: 'copy'
    
    input:
    path samplesheet
    path master_ped
    
    output:
    path "samplesheet.valid.csv", emit: csv
    path "batch.ped", emit: batch_ped
    path "master.ped", emit: master_ped
    path "master.ped_backup_*", emit: backup_ped, optional: true
    path "versions.yml", emit: versions
    
    script:
    def has_master = master_ped.name != 'NO_FILE'
    """
    python3 << 'PYTHON_SCRIPT'
import csv
import sys
import os
import shutil
from pathlib import Path
from datetime import datetime

# Validate CSV with standard PED order:
# family_id,sample_id,paternal_id,maternal_id,sex,phenotype,fastq_r1,fastq_r2
def validate_samplesheet(samplesheet_path):
    required_cols = ['family_id', 'sample_id', 'fastq_r1', 'fastq_r2']
    optional_cols = ['paternal_id', 'maternal_id', 'sex', 'phenotype']
    
    samples = []
    errors = []
    
    with open(samplesheet_path, 'r') as f:
        reader = csv.DictReader(f)
        headers = reader.fieldnames
        
        print(f"CSV columns found: {headers}")
        print(f"Expected order: family_id,sample_id,paternal_id,maternal_id,sex,phenotype,fastq_r1,fastq_r2")
        
        for col in required_cols:
            if col not in headers:
                errors.append(f"Missing required column: {col}")
        
        if errors:
            print("\\n".join(errors), file=sys.stderr)
            sys.exit(1)
        
        for i, row in enumerate(reader, start=2):
            sample = {}
            
            sample['family_id'] = row['family_id'].strip()
            sample['sample_id'] = row['sample_id'].strip()
            sample['fastq_r1'] = row['fastq_r1'].strip()
            sample['fastq_r2'] = row['fastq_r2'].strip()
            
            if not sample['family_id']:
                errors.append(f"Row {i}: Empty family_id")
            if not sample['sample_id']:
                errors.append(f"Row {i}: Empty sample_id")
            if not sample['fastq_r1'] or not os.path.exists(sample['fastq_r1']):
                errors.append(f"Row {i}: fastq_r1 not found: {sample['fastq_r1']}")
            if not sample['fastq_r2'] or not os.path.exists(sample['fastq_r2']):
                errors.append(f"Row {i}: fastq_r2 not found: {sample['fastq_r2']}")
            
            sample['paternal_id'] = row.get('paternal_id', '0').strip() or '0'
            sample['maternal_id'] = row.get('maternal_id', '0').strip() or '0'
            
            sex = row.get('sex', '0').strip().lower()
            if sex in ['male', 'm', '1']:
                sample['sex'] = '1'
            elif sex in ['female', 'f', '2']:
                sample['sex'] = '2'
            else:
                sample['sex'] = '0'
            
            pheno = row.get('phenotype', '0').strip().lower()
            if pheno in ['affected', 'a', '2']:
                sample['phenotype'] = '2'
            elif pheno in ['unaffected', 'u', '1']:
                sample['phenotype'] = '1'
            else:
                sample['phenotype'] = '0'
            
            samples.append(sample)
    
    if errors:
        print("Validation errors:", file=sys.stderr)
        print("\\n".join(errors), file=sys.stderr)
        sys.exit(1)
    
    return samples

# Write validated CSV in standard PED order
def write_valid_csv(samples, output_path):
    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=[
            'family_id', 'sample_id', 'paternal_id', 'maternal_id', 
            'sex', 'phenotype', 'fastq_r1', 'fastq_r2'
        ])
        writer.writeheader()
        writer.writerows(samples)

# Write PED file (6 columns: FID, IID, PAT, MAT, SEX, PHENO)
def write_ped(samples, output_path):
    with open(output_path, 'w') as f:
        for s in samples:
            line = f"{s['family_id']}\\t{s['sample_id']}\\t{s['paternal_id']}\\t{s['maternal_id']}\\t{s['sex']}\\t{s['phenotype']}\\n"
            f.write(line)

# Read existing PED file
def read_existing_ped(ped_path):
    samples = {}
    if os.path.exists(ped_path) and os.path.getsize(ped_path) > 0:
        with open(ped_path, 'r') as f:
            for line in f:
                if line.strip():
                    parts = line.strip().split('\\t')
                    if len(parts) >= 6:
                        sample_id = parts[1]
                        samples[sample_id] = {
                            'family_id': parts[0],
                            'sample_id': parts[1],
                            'paternal_id': parts[2],
                            'maternal_id': parts[3],
                            'sex': parts[4],
                            'phenotype': parts[5]
                        }
    return samples

# Create timestamped backup of existing PED file
def backup_ped(ped_path):
    if os.path.exists(ped_path) and os.path.getsize(ped_path) > 0:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        backup_path = f"master.ped_backup_{timestamp}"
        shutil.copy(ped_path, backup_path)
        print(f"Created backup: {backup_path}")
        return backup_path
    return None

# Merge new samples into existing PED, updating if exists
def merge_ped(existing_samples, new_samples):
    merged = existing_samples.copy()
    for s in new_samples:
        merged[s['sample_id']] = s
    return list(merged.values())

# Main execution
print("=== Input Validation ===")

samples = validate_samplesheet('${samplesheet}')
print(f"Validated {len(samples)} samples")

write_valid_csv(samples, 'samplesheet.valid.csv')
print("Written: samplesheet.valid.csv")

write_ped(samples, 'batch.ped')
print("Written: batch.ped")

has_master = ${has_master ? 'True' : 'False'}
if has_master:
    backup_path = backup_ped('${master_ped}')
    
    existing = read_existing_ped('${master_ped}')
    print(f"Read {len(existing)} samples from existing master.ped")
    merged = merge_ped(existing, samples)
    write_ped(merged, 'master.ped')
    print(f"Written merged master.ped with {len(merged)} samples")
else:
    write_ped(samples, 'master.ped')
    print(f"Written new master.ped with {len(samples)} samples")

print("=== Validation Complete ===")
PYTHON_SCRIPT
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}

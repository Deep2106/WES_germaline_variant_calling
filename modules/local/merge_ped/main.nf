/*
========================================================================================
    MERGE PED MODULE
========================================================================================
    Merges new sample PED entries with existing master PED file.
    Used for incremental joint calling to maintain a complete pedigree.

    Changes vs original:
    - Deduplication is now explicit: last-write-wins per (family_id, sample_id) key.
      This handles the case where INPUT_CHECK may emit one PED line per lane
      before dedup occurs upstream, or if a sample is resubmitted across batches.

    Input:
        new_ped      - PED file from current batch (already deduped by INPUT_CHECK)
        existing_ped - Existing master PED file (optional, may be NO_FILE)
        output_path  - Persistent path to copy master.ped (optional)

    Output:
        master.ped   - Combined, deduplicated PED file
        versions.yml

    PED columns (tab-separated): FamilyID  IndividualID  PaternalID  MaternalID  Sex  Phenotype
    Sex: 1=male, 2=female, 0=unknown | Phenotype: 1=unaffected, 2=affected, 0=unknown
----------------------------------------------------------------------------------------
*/

process MERGE_PED {
    tag "merge_ped"
    label 'process_single'
    container "${params.containers?.python ?: 'python:3.9-slim'}"

    publishDir "${params.outdir}/pedigree", mode: params.publish_dir_mode ?: 'copy', failOnError: false

    input:
    path new_ped
    path existing_ped
    val  output_path

    output:
    path "master.ped",   emit: merged_ped
    path "versions.yml", emit: versions

    script:
    def has_existing = existing_ped.name != 'NO_FILE' && existing_ped.name != 'NO_PED'
    """
    python3 << 'PYTHON_SCRIPT'
import sys

# ---------------------------------------------------------------------------
# Read PED into ordered dict keyed by (family_id, sample_id).
# Last entry wins for duplicates (new batch overrides existing).
# Comments (lines starting with #) are preserved in insertion order.
# ---------------------------------------------------------------------------
def read_ped(path):
    entries  = {}   # (fam, ind) -> line
    comments = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\\n")
            if not line.strip():
                continue
            if line.startswith('#'):
                comments.append(line)
                continue
            parts = line.split('\\t')
            if len(parts) < 6:
                print(f"WARNING: Skipping malformed line: {line!r}", file=sys.stderr)
                continue
            # Validate sex and phenotype
            if parts[4] not in ('0','1','2'):
                print(f"WARNING: Invalid sex value '{parts[4]}' for {parts[1]}", file=sys.stderr)
            if parts[5] not in ('0','1','2'):
                print(f"WARNING: Invalid phenotype '{parts[5]}' for {parts[1]}", file=sys.stderr)
            key = (parts[0], parts[1])
            entries[key] = '\\t'.join(parts[:6])
    return comments, entries

def write_ped(comments, entries, path):
    with open(path, 'w') as f:
        for c in comments:
            f.write(c + '\\n')
        for line in entries.values():
            f.write(line + '\\n')

print("=== MERGE PED ===")

comments = []
merged   = {}

has_existing = ${has_existing ? 'True' : 'False'}

# Process existing master first (so new batch can override)
if has_existing:
    c, e = read_ped('${existing_ped}')
    comments.extend(c)
    merged.update(e)
    print(f"Existing master.ped: {len(e)} sample(s)")
else:
    print("No existing master.ped (first run)")

# Process new batch (overrides duplicates)
c_new, e_new = read_ped('${new_ped}')
# Only add comments not already present
for c in c_new:
    if c not in comments:
        comments.append(c)
overlap = len(set(e_new) & set(merged))
merged.update(e_new)
print(f"New batch: {len(e_new)} sample(s), {overlap} updated existing entry/entries")

write_ped(comments, merged, 'master.ped')

total_samples  = len(merged)
total_families = len({k[0] for k in merged})
print(f"master.ped: {total_samples} sample(s) across {total_families} family/families")

# Copy to persistent path if provided
output_path = '${output_path}'
if output_path and output_path != 'null':
    import os, shutil
    os.makedirs(os.path.dirname(output_path) or '.', exist_ok=True)
    shutil.copy('master.ped', output_path)
    print(f"Copied to persistent path: {output_path}")

print("=== MERGE PED COMPLETE ===")

# Write versions
with open('versions.yml', 'w') as f:
    import subprocess
    py_ver = subprocess.check_output(['python3', '--version']).decode().strip().replace('Python ', '')
    f.write(f'"${task.process}":\\n')
    f.write(f'    python: {py_ver}\\n')
    f.write(f'    total_samples: {total_samples}\\n')
    f.write(f'    total_families: {total_families}\\n')

PYTHON_SCRIPT
    """

    stub:
    """
    echo -e "FAM001\\tSAMPLE001\\t0\\t0\\t1\\t2" > master.ped

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}

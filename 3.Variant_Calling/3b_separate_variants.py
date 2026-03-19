#!/usr/bin/env python3
"""
Step 3b: Separate SNVs and Indels from Mutect2 VCFs

Splits Mutect2 VCF files into separate SNV and indel files for downstream filtering.
Submits one SGE job per stringency level.

SNVs are extracted from all three stringency levels (MS1, MS2, Duplex).
Indels are extracted from Duplex only.

Input:
    - Mutect2 VCFs in MS1/With_rg/Mutect_output, MS2/With_rg/Mutect_output,
      DUPLEX/With_rg/Mutect_output

Output:
    - {base_dir}/Candidates/SNVs-MS1/
    - {base_dir}/Candidates/SNVs-MS2/
    - {base_dir}/Candidates/SNVs-duplex/
    - {base_dir}/Candidates/Indels/

Requirements:
    - SGE cluster with qsub
"""

import os
import subprocess
import argparse
from pathlib import Path


def find_vcf_files(directory):
    """Find all VCF files in a directory."""
    vcf_files = []
    for filename in os.listdir(directory):
        if filename.endswith(".vcf") or filename.endswith(".vcf.gz"):
            vcf_files.append(os.path.join(directory, filename))
    return sorted(vcf_files)


def generate_split_script(vcf_file, snv_output, indel_output, process_indels):
    """
    Generate inline Python to split one VCF into SNVs and indels.
    Returns a string of Python code for embedding in the SGE job script.
    """
    indel_arg = f'"{indel_output}"' if process_indels else 'None'

    script = f"""
import gzip, sys

input_vcf = "{vcf_file}"
snv_output = "{snv_output}"
indel_output = {indel_arg}

snv_count = 0
indel_count = 0

opener = gzip.open if input_vcf.endswith('.gz') else open
mode = 'rt' if input_vcf.endswith('.gz') else 'r'

with opener(input_vcf, mode) as infile, open(snv_output, 'w') as snv_file:
    indel_file = open(indel_output, 'w') if indel_output else None
    for line in infile:
        if line.startswith('#'):
            snv_file.write(line)
            if indel_file:
                indel_file.write(line)
            continue
        cols = line.strip().split('\\t')
        ref, alt = cols[3], cols[4]
        if len(ref) == 1 and len(alt) == 1:
            snv_file.write(line)
            snv_count += 1
        else:
            if indel_file:
                indel_file.write(line)
            indel_count += 1
    if indel_file:
        indel_file.close()

print(f"  SNVs written: {{snv_count}}")
if indel_output:
    print(f"  Indels written: {{indel_count}}")
else:
    print(f"  Indels skipped: {{indel_count}}")
"""
    return script


def create_job_script(label, vcf_files, output_dir, job_script_dir, sge_config, process_indels):
    """
    Create an SGE job script to split all VCFs for one stringency level.

    Returns:
        Path to created job script.
    """
    snv_dir = os.path.join(output_dir, f"SNVs-{label}")
    indel_dir = os.path.join(output_dir, "Indels") if process_indels else None

    log_file = os.path.join(job_script_dir, f"sep05b_{label}_separate.log")
    err_file = os.path.join(job_script_dir, f"sep05b_{label}_separate.err")
    job_script_path = os.path.join(job_script_dir, f"sep05b_{label}_separate.sh")

    python_blocks = []
    for vcf_file in vcf_files:
        sample_name = os.path.basename(vcf_file).replace('.vcf.gz', '').replace('.vcf', '')
        snv_output = os.path.join(snv_dir, f"{sample_name}.vcf")
        indel_output = os.path.join(indel_dir, f"{sample_name}.vcf") if process_indels else None

        block = f'print("Processing: {sample_name}")\n'
        block += generate_split_script(vcf_file, snv_output, indel_output, process_indels)
        python_blocks.append(block)

    all_python = "\n".join(python_blocks)

    job_script = f"""#!/bin/bash -l
#$ -N sep05b_{label}
#$ -l h_rt={sge_config.get('wall_time', '01:00:00')}
#$ -l h_vmem={sge_config.get('memory', '8G')}
#$ -pe smp {sge_config.get('cores', 1)}
#$ -o {log_file}
#$ -e {err_file}
#$ -cwd

echo "============================================================"
echo "Step 3b: Separate SNVs and Indels - {label.upper()}"
echo "============================================================"
echo "Started: $(date)"
echo ""

mkdir -p "{snv_dir}"
"""
    if process_indels:
        job_script += f'mkdir -p "{indel_dir}"\n'

    job_script += f"""
echo "Output SNVs → {snv_dir}"
"""
    if process_indels:
        job_script += f'echo "Output Indels → {indel_dir}"\n'

    job_script += f"""
echo ""

python3 - << 'PYEOF'
{all_python}
PYEOF

echo ""
echo "============================================================"
echo "Step 3b {label.upper()} complete: $(date)"
echo "============================================================"
"""

    with open(job_script_path, 'w') as f:
        f.write(job_script)
    os.chmod(job_script_path, 0o755)

    return job_script_path


def main():
    parser = argparse.ArgumentParser(
        description='Separate SNVs and indels from Mutect2 VCF files (SGE submission)',
        epilog='Example: python 3b_separate_variants.py --base-dir /path/to/creak_output'
    )
    parser.add_argument(
        '--base-dir',
        required=True,
        help='Base directory containing MS1/, MS2/, DUPLEX/ subdirectories'
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Show what would be processed without submitting jobs'
    )

    args = parser.parse_args()

    sge_config = {
        'memory': '8G',
        'cores': 1,
        'wall_time': '01:00:00'
    }

    base_dir = os.path.abspath(args.base_dir)
    if not os.path.isdir(base_dir):
        print(f"ERROR: Base directory not found: {base_dir}")
        return

    stringencies = {
        'MS1':    (os.path.join(base_dir, 'MS1',    'With_rg', 'Mutect_output'), False),
        'MS2':    (os.path.join(base_dir, 'MS2',    'With_rg', 'Mutect_output'), False),
        'duplex': (os.path.join(base_dir, 'DUPLEX', 'With_rg', 'Mutect_output'), True),
    }

    available = {}
    for label, (path, process_indels) in stringencies.items():
        if os.path.isdir(path):
            vcfs = find_vcf_files(path)
            if vcfs:
                available[label] = (path, process_indels, vcfs)
            else:
                print(f"Warning: No VCF files found in {path}")
        else:
            print(f"Warning: Directory not found: {path}")

    if not available:
        print("ERROR: No Mutect_output directories with VCF files found.")
        return

    output_dir = os.path.join(base_dir, 'Candidates')
    job_script_dir = os.path.join(base_dir, 'job_scripts', '05b_separate')

    print(f"\n{'='*60}")
    print(f"Step 3b: Separate SNVs and Indels")
    print(f"{'='*60}")
    print(f"Base directory:  {base_dir}")
    print(f"Output:          {output_dir}/")
    print(f"Job scripts:     {job_script_dir}/")
    print(f"Stringencies:    {', '.join(available.keys())}")
    print()

    if args.dry_run:
        print("DRY RUN - no jobs submitted\n")
        for label, (path, process_indels, vcfs) in available.items():
            print(f"  [{label}] {len(vcfs)} VCF(s) found")
            print(f"    Input:  {path}")
            print(f"    SNVs → Candidates/SNVs-{label}/")
            if process_indels:
                print(f"    Indels → Candidates/Indels/")
        return

    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(job_script_dir, exist_ok=True)

    submitted_jobs = []

    for label, (path, process_indels, vcfs) in available.items():
        print(f"Creating job for {label.upper()} ({len(vcfs)} VCF files)...")

        job_script_path = create_job_script(
            label=label,
            vcf_files=vcfs,
            output_dir=output_dir,
            job_script_dir=job_script_dir,
            sge_config=sge_config,
            process_indels=process_indels
        )

        result = subprocess.run(['qsub', job_script_path], capture_output=True, text=True)
        if result.returncode == 0:
            job_id = result.stdout.strip().split()[2]
            submitted_jobs.append(job_id)
            print(f"  ✓ Submitted: {os.path.basename(job_script_path)}")
            print(f"    Job ID: {job_id}")
        else:
            print(f"  ✗ Submission failed: {result.stderr}")

    print(f"\n{'='*60}")
    print(f"Submitted {len(submitted_jobs)} job(s): {', '.join(submitted_jobs)}")
    print(f"Monitor with: qstat -u $USER")
    print(f"{'='*60}\n")

    print("When complete, check outputs:")
    print(f"  ls -lh {output_dir}/SNVs-MS1/")
    print(f"  ls -lh {output_dir}/SNVs-MS2/")
    print(f"  ls -lh {output_dir}/SNVs-duplex/")
    print(f"  ls -lh {output_dir}/Indels/")


if __name__ == "__main__":
    main()

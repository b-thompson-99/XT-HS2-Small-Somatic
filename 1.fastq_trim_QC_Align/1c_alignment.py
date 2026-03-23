#!/usr/bin/env python3
"""
Step 1c: BWA-MEM Alignment with MBC Tag Preservation

Generates and submits SGE jobs to align paired-end FASTQ files to a reference
genome using BWA-MEM. The -C flag preserves molecular barcode (MBC) tags in
the BAM for downstream consensus collapse.

Input:
    - Directory containing sample subdirectories, each with QC'd FASTQ files
      ending in *_R1_fastp.fastq.gz and *_R2_fastp.fastq.gz

Output:
    - Sorted BAMs with MBC tags: {sample}_sorted.bam

Requirements:
    - BWA (v0.7.17 or later)
    - samtools (v1.10 or later)
    - Reference genome with BWA index
    - SGE cluster with qsub
    - config_1c_alignment.yaml
"""

import os
import subprocess
import time
import yaml
import argparse
from pathlib import Path


def load_config(config_path):
    """Load configuration from YAML file."""
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)


def find_fastq_pairs(sample_dir, r1_pattern, r2_pattern):
    """
    Find paired FASTQ files in a sample directory.

    Returns:
        Tuple of (R1_path, R2_path) or (None, None) if not found.
    """
    r1_candidates = [f for f in os.listdir(sample_dir) if f.endswith(r1_pattern)]
    r2_candidates = [f for f in os.listdir(sample_dir) if f.endswith(r2_pattern)]

    if len(r1_candidates) != 1 or len(r2_candidates) != 1:
        print(f"Warning: Expected one R1 and one R2 file in {sample_dir}")
        print(f"  R1 candidates: {r1_candidates}")
        print(f"  R2 candidates: {r2_candidates}")
        return None, None

    r1_path = os.path.join(sample_dir, r1_candidates[0])
    r2_path = os.path.join(sample_dir, r2_candidates[0])

    return r1_path, r2_path


def verify_dependencies(config):
    """
    Verify required tools and reference genome exist.

    Raises:
        FileNotFoundError: If dependencies are missing.
    """
    bwa_path = config['tools']['bwa_path']
    samtools_path = config['tools']['samtools_path']
    reference = config['reference']['genome_fasta']

    if not os.path.exists(bwa_path):
        raise FileNotFoundError(f"BWA not found: {bwa_path}")

    if not os.path.exists(samtools_path):
        raise FileNotFoundError(f"samtools not found: {samtools_path}")

    if not os.path.exists(reference):
        raise FileNotFoundError(f"Reference genome not found: {reference}")

    index_extensions = ['.amb', '.ann', '.bwt', '.pac', '.sa']
    missing_indices = [ext for ext in index_extensions if not os.path.exists(reference + ext)]

    if missing_indices:
        raise FileNotFoundError(
            f"BWA index incomplete for {reference}\n"
            f"Missing extensions: {missing_indices}\n"
            f"Run: bwa index {reference}"
        )

    print(f"✓ BWA found: {bwa_path}")
    print(f"✓ samtools found: {samtools_path}")
    print(f"✓ Reference genome found: {reference}")
    print(f"✓ BWA index complete")


def create_qsub_script(sample_name, config, fastq1, fastq2, output_dir, script_dir):
    """
    Create SGE job script for BWA-MEM alignment.

    Returns:
        Path to created job script.
    """
    output_bam = os.path.join(output_dir, f"{sample_name}_sorted.bam")

    job_script_template = """#!/bin/bash
#$ -N align_{sample}
#$ -l h_vmem={memory}
#$ -pe smp {cores}
#$ -l h_rt={wall_time}
#$ -cwd
#$ -o {output_dir}/{sample}_align.stdout
#$ -e {output_dir}/{sample}_align.stderr

# Verify tools are available
which bwa
which samtools
bwa 2>&1 | head -n 3
samtools --version | head -n 1

# Input files
FQ1={fq1}
FQ2={fq2}

# Output file
SORTED_BAM={output_bam}

# Reference genome
REF={reference}

echo "Starting alignment for {sample}..."
echo "  Input R1: $(basename $FQ1)"
echo "  Input R2: $(basename $FQ2)"
echo "  Reference: $(basename $REF)"

# Run BWA-MEM with -C flag to preserve MBC tags
# Then pipe to samtools for BAM conversion and sorting
{bwa_path} mem -C -t {cores} \\
  -R '@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINA' \\
  $REF $FQ1 $FQ2 | \\
  {samtools_path} view -b - | \\
  {samtools_path} sort -@ {sort_threads} -o $SORTED_BAM

# Verify output
if [ -f "$SORTED_BAM" ]; then
    SIZE=$(stat -f%z "$SORTED_BAM" 2>/dev/null || stat -c%s "$SORTED_BAM" 2>/dev/null)
    echo "✓ Alignment completed successfully"
    echo "  Output BAM: $(basename $SORTED_BAM)"
    echo "  File size: $SIZE bytes"
    echo "  Read count: $({samtools_path} view -c $SORTED_BAM)"
else
    echo "✗ ERROR: Output BAM not created"
    exit 1
fi
"""

    sort_threads = max(1, config['sge']['cores'] - 1)

    job_script = job_script_template.format(
        sample=sample_name,
        memory=config['sge']['memory'],
        cores=config['sge']['cores'],
        wall_time=config['sge']['wall_time'],
        output_dir=output_dir,
        bwa_path=config['tools']['bwa_path'],
        samtools_path=config['tools']['samtools_path'],
        reference=config['reference']['genome_fasta'],
        fq1=fastq1,
        fq2=fastq2,
        output_bam=output_bam,
        sort_threads=sort_threads
    )

    job_script_path = os.path.join(script_dir, f"{sample_name}_align.sh")
    with open(job_script_path, "w") as f:
        f.write(job_script)

    return job_script_path


def submit_jobs_in_batches(job_scripts, batch_size, pause_time):
    """Submit jobs to SGE in batches."""
    for i in range(0, len(job_scripts), batch_size):
        batch = job_scripts[i:i + batch_size]

        print(f"\n=== Submitting batch {i // batch_size + 1} ({len(batch)} jobs) ===")

        for job_script in batch:
            result = subprocess.run(
                ["qsub", job_script],
                capture_output=True,
                text=True
            )
            if result.returncode == 0:
                print(f"  ✓ Submitted: {os.path.basename(job_script)}")
                print(f"    {result.stdout.strip()}")
            else:
                print(f"  ✗ Failed: {os.path.basename(job_script)}")
                print(f"    Error: {result.stderr.strip()}")

        if i + batch_size < len(job_scripts):
            print(f"  Pausing {pause_time}s before next batch...")
            time.sleep(pause_time)


def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(
        description='Generate and submit BWA-MEM alignment jobs with MBC tag preservation'
    )
    parser.add_argument(
        '--config',
        default='config/03_alignment_config.yaml',
        help='Path to configuration file (default: config/03_alignment_config.yaml)'
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Generate job scripts but do not submit them'
    )

    args = parser.parse_args()

    config = load_config(args.config)

    print("Checking dependencies...")
    verify_dependencies(config)
    print()

    input_dir = config['paths']['input_fastq_dir']
    script_dir = config['paths']['job_scripts_dir']

    os.makedirs(script_dir, exist_ok=True)

    job_scripts = []

    print(f"Scanning for samples in: {input_dir}\n")

    for sample in os.listdir(input_dir):
        sample_dir = os.path.join(input_dir, sample)

        if not os.path.isdir(sample_dir):
            continue

        print(f"Processing sample: {sample}")

        fastq1, fastq2 = find_fastq_pairs(
            sample_dir,
            config['fastq']['r1_suffix'],
            config['fastq']['r2_suffix']
        )

        if fastq1 is None or fastq2 is None:
            print(f"  ⚠ Skipping {sample} - FASTQ files not found\n")
            continue

        print(f"  R1: {os.path.basename(fastq1)}")
        print(f"  R2: {os.path.basename(fastq2)}")

        job_script_path = create_qsub_script(
            sample_name=sample,
            config=config,
            fastq1=fastq1,
            fastq2=fastq2,
            output_dir=sample_dir,
            script_dir=script_dir
        )

        job_scripts.append(job_script_path)
        print(f"  ✓ Job script created: {os.path.basename(job_script_path)}\n")

    print(f"{'='*60}")
    print(f"Generated {len(job_scripts)} job scripts")
    print(f"{'='*60}\n")

    if len(job_scripts) == 0:
        print("No jobs to submit. Exiting.")
        return

    if args.dry_run:
        print("DRY RUN: Jobs not submitted.")
        print("Job scripts saved to:", script_dir)
    else:
        submit_jobs_in_batches(
            job_scripts,
            batch_size=config['submission']['batch_size'],
            pause_time=config['submission']['pause_between_batches']
        )
        print("\n✓ All batches submitted successfully")


if __name__ == "__main__":
    main()

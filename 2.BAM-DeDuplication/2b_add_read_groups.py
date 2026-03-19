#!/usr/bin/env python3
"""
Step 2b: Add Read Groups to Consensus BAMs

Generates and submits SGE jobs to add read group (RG) tags to consensus BAMs
from step 2a. RG tags are required by Mutect2 for variant calling.

Processes all three stringency levels (MS1, MS2, DUPLEX) automatically.

Input:
    - Consensus BAMs from step 2a in MS1/, MS2/, DUPLEX/ subdirectories

Output:
    - RG-tagged BAMs in With_rg/ subdirectories
    - Indexed BAM files (.bam.bai)

Requirements:
    - Picard
    - samtools (v1.10 or later)
    - Java 17+
    - SGE cluster with qsub
    - config_2b_add_read_groups.yaml
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


def find_bam_files(directory):
    """Find all BAM files in a directory, excluding already-processed with_rg BAMs."""
    bam_files = []
    for filename in sorted(os.listdir(directory)):
        if filename.endswith('.bam') and not filename.endswith('.with_rg.bam'):
            bam_files.append(os.path.join(directory, filename))
    return bam_files


def extract_sample_name(bam_filename):
    """
    Extract sample name from BAM filename.

    Examples:
        SAMPLE_001_hybrid_MS1.bam → SAMPLE_001
        SAMPLE_002_hybrid_MS2.bam → SAMPLE_002
        SAMPLE_003_duplex.bam     → SAMPLE_003
    """
    basename = os.path.basename(bam_filename).replace('.bam', '')
    basename = basename.replace('_hybrid_MS1', '')
    basename = basename.replace('_hybrid_MS2', '')
    basename = basename.replace('_duplex', '')
    return basename


def verify_dependencies(config):
    """
    Verify required tools exist.

    Raises:
        FileNotFoundError: If dependencies are missing.
    """
    picard_path = config['tools']['picard']
    samtools_path = config['tools']['samtools']

    if not os.path.exists(picard_path):
        raise FileNotFoundError(f"Picard not found: {picard_path}")
    if not os.path.exists(samtools_path):
        raise FileNotFoundError(f"samtools not found: {samtools_path}")

    print(f"✓ Picard found: {picard_path}")
    print(f"✓ samtools found: {samtools_path}")


def create_qsub_script(sample_name, stringency_label, input_bam, output_bam, config, script_dir):
    """
    Create SGE job script for Picard AddOrReplaceReadGroups + samtools index.

    Returns:
        Path to created job script.
    """
    rg = config['read_group_params']
    picard = config['tools']['picard']
    samtools = config['tools']['samtools']
    java_home = config['tools']['java_home']
    output_dir = os.path.dirname(output_bam)

    job_script_template = """#!/bin/bash
#$ -N add_rg_{sample}_{stringency}
#$ -l h_vmem={memory}
#$ -pe smp {cores}
#$ -l h_rt={wall_time}
#$ -cwd
#$ -o {output_dir}/{sample}_add_rg_{stringency}.stdout
#$ -e {output_dir}/{sample}_add_rg_{stringency}.stderr

# Set Java environment (required for Picard)
export JAVA_HOME={java_home}
export PATH={java_home}/bin:$PATH

INPUT_BAM={input_bam}
OUTPUT_BAM={output_bam}

echo "Adding read groups for {sample} ({stringency})..."
echo "  Input:  $(basename $INPUT_BAM)"
echo "  Output: $(basename $OUTPUT_BAM)"

# Add read groups with Picard
{picard} AddOrReplaceReadGroups \\
    I=$INPUT_BAM \\
    O=$OUTPUT_BAM \\
    RGID={rgid} \\
    RGLB={rglb} \\
    RGPL={rgpl} \\
    RGPU={rgpu} \\
    RGSM={sample}

if [ $? -ne 0 ]; then
    echo "✗ ERROR: Picard failed for {sample}"
    exit 1
fi

# Index BAM
echo "  Indexing BAM..."
{samtools} index $OUTPUT_BAM

if [ $? -ne 0 ]; then
    echo "✗ ERROR: Indexing failed for {sample}"
    exit 1
fi

echo "✓ Read groups added and indexed successfully"
echo "  Output: $(basename $OUTPUT_BAM)"
"""

    job_script = job_script_template.format(
        sample=sample_name,
        stringency=stringency_label,
        memory=config['sge']['memory'],
        cores=config['sge']['cores'],
        wall_time=config['sge']['wall_time'],
        output_dir=output_dir,
        java_home=java_home,
        input_bam=input_bam,
        output_bam=output_bam,
        picard=picard,
        samtools=samtools,
        rgid=rg.get('rgid', '1'),
        rglb=rg.get('rglb', 'lib1'),
        rgpl=rg.get('rgpl', 'ILLUMINA'),
        rgpu=rg.get('rgpu', 'unit1')
    )

    job_script_path = os.path.join(script_dir, f"{sample_name}_add_rg_{stringency_label}.sh")
    with open(job_script_path, 'w') as f:
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
        description='Generate and submit SGE jobs to add read groups to CReaK consensus BAMs'
    )
    parser.add_argument(
        '--config',
        required=True,
        help='Path to YAML configuration file'
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

    base_dir = config['paths']['creak_output_dir']
    script_dir = config['paths']['job_scripts_dir']

    os.makedirs(script_dir, exist_ok=True)

    stringencies = ['MS1', 'MS2', 'DUPLEX']
    all_job_scripts = []

    for stringency in stringencies:
        stringency_dir = os.path.join(base_dir, stringency)
        output_dir = os.path.join(stringency_dir, 'With_rg')

        print(f"\n{'='*60}")
        print(f"Processing {stringency} BAMs")
        print(f"{'='*60}\n")

        if not os.path.isdir(stringency_dir):
            print(f"  ⚠ Directory not found: {stringency_dir} - skipping\n")
            continue

        bam_files = find_bam_files(stringency_dir)

        if not bam_files:
            print(f"  ⚠ No BAM files found in {stringency_dir} - skipping\n")
            continue

        print(f"Found {len(bam_files)} BAM file(s)")
        print(f"Output directory: {output_dir}\n")

        os.makedirs(output_dir, exist_ok=True)

        for input_bam in bam_files:
            basename = os.path.basename(input_bam)
            sample_name = extract_sample_name(input_bam)
            output_basename = basename.replace('.bam', '.with_rg.bam')
            output_bam = os.path.join(output_dir, output_basename)

            print(f"  Sample:  {sample_name}")
            print(f"  Input:   {basename}")
            print(f"  Output:  {output_basename}")

            job_script_path = create_qsub_script(
                sample_name=sample_name,
                stringency_label=stringency,
                input_bam=input_bam,
                output_bam=output_bam,
                config=config,
                script_dir=script_dir
            )

            all_job_scripts.append(job_script_path)
            print(f"  ✓ Job script created: {os.path.basename(job_script_path)}\n")

    print(f"\n{'='*60}")
    print(f"Generated {len(all_job_scripts)} job scripts total")
    print(f"{'='*60}\n")

    if not all_job_scripts:
        print("No jobs to submit. Exiting.")
        return

    if args.dry_run:
        print("DRY RUN: Jobs not submitted.")
        print(f"Job scripts saved to: {script_dir}")
    else:
        submit_jobs_in_batches(
            all_job_scripts,
            batch_size=config['submission']['batch_size'],
            pause_time=config['submission']['pause_between_batches']
        )
        print("\n✓ All batches submitted successfully")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Step 1b: FASTQ QC with fastp

Generates and submits SGE jobs to run fastp QC on paired-end FASTQ files.

Input:
    - Directory containing sample subdirectories, each with paired FASTQ files
      ending in *_trimmed_R1.fastq.gz and *_trimmed_R2.fastq.gz

Output:
    - QC'd FASTQs: {sample}_trimmed_R1/R2_fastp.fastq.gz
    - QC reports: {sample}_fastp_report.html / .json

Requirements:
    - fastp
    - SGE cluster with qsub
    - config_1b_fastq_qc.yaml
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


def create_qsub_script(sample_name, config, fastq1, fastq2, output_dir, script_dir):
    """
    Create SGE job script for fastp QC.

    Returns:
        Path to created job script.
    """
    output_r1 = os.path.join(output_dir, f"{sample_name}_trimmed_R1_fastp.fastq.gz")
    output_r2 = os.path.join(output_dir, f"{sample_name}_trimmed_R2_fastp.fastq.gz")
    html_report = os.path.join(output_dir, f"{sample_name}_fastp_report.html")
    json_report = os.path.join(output_dir, f"{sample_name}_fastp_report.json")

    job_script_template = """#!/bin/bash
#$ -N fastp_qc_{sample}
#$ -l h_vmem={memory}
#$ -pe smp {cores}
#$ -l h_rt={wall_time}
#$ -cwd
#$ -o {output_dir}/{sample}_fastp.stdout
#$ -e {output_dir}/{sample}_fastp.stderr

# Input files
R1={fq1}
R2={fq2}

# Output files
OUTPUT_R1={output_r1}
OUTPUT_R2={output_r2}
HTML_REPORT={html_report}
JSON_REPORT={json_report}

# Run fastp
{fastp_path} -i $R1 -I $R2 \\
      -o $OUTPUT_R1 -O $OUTPUT_R2 \\
      -h $HTML_REPORT -j $JSON_REPORT \\
      --thread {cores}

echo "fastp QC completed for {sample}"
"""

    job_script = job_script_template.format(
        sample=sample_name,
        memory=config['sge']['memory'],
        cores=config['sge']['cores'],
        wall_time=config['sge']['wall_time'],
        output_dir=output_dir,
        fastp_path=config['tools']['fastp_path'],
        fq1=fastq1,
        fq2=fastq2,
        output_r1=output_r1,
        output_r2=output_r2,
        html_report=html_report,
        json_report=json_report
    )

    job_script_path = os.path.join(script_dir, f"{sample_name}_fastp_qc.sh")
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
        description='Generate and submit SGE jobs for fastp QC'
    )
    parser.add_argument(
        '--config',
        default='config/02_fastq_qc_config.yaml',
        help='Path to configuration file (default: config/02_fastq_qc_config.yaml)'
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Generate job scripts but do not submit them'
    )

    args = parser.parse_args()

    config = load_config(args.config)

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

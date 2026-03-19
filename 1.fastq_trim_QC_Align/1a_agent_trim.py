#!/usr/bin/env python3
"""
Step 1a: MBC Trimming with AGeNT

Generates and submits SGE jobs to run Agilent AGeNT trimming on paired-end FASTQ files.
Extracts molecular barcodes (MBCs/UMIs) and prepares reads for duplex consensus calling.

Input:
    - Directory containing sample subdirectories, each with paired FASTQ files
      ending in *_R1.fastq.gz and *_R2.fastq.gz

Output:
    - Trimmed FASTQs with MBC tags: {sample}_trimmed_R1/R2.fastq.gz
    - AGeNT trim reports and logs

Requirements:
    - AGeNT (v3.1.2 or later)
    - Java 17+
    - SGE cluster with qsub
    - config_1a_agent_trim.yaml
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
    Verify required tools exist.

    Raises:
        FileNotFoundError: If dependencies are missing.
    """
    agent_path = config['tools']['agent_path']
    java_home = config['tools']['java_home']

    if not os.path.exists(agent_path):
        raise FileNotFoundError(f"AGeNT tool not found: {agent_path}")

    if not os.path.exists(java_home):
        raise FileNotFoundError(f"Java installation not found: {java_home}")

    print(f"✓ AGeNT found: {agent_path}")
    print(f"✓ Java found: {java_home}")


def create_qsub_script(sample_name, config, fastq1, fastq2, output_dir, script_dir):
    """
    Create SGE job script for AGeNT trimming.

    Returns:
        Path to created job script.
    """
    output_prefix = os.path.join(output_dir, f"{sample_name}_trimmed")

    java_home = config['tools']['java_home']
    java_env = f"""export JAVA_HOME={java_home}
export PATH={java_home}/bin:$PATH"""

    job_script_template = """#!/bin/bash
#$ -N agent_trim_{sample}
#$ -l h_vmem={memory}
#$ -pe smp {cores}
#$ -l h_rt={wall_time}
#$ -cwd
#$ -o {output_dir}/{sample}_trim.stdout
#$ -e {output_dir}/{sample}_trim.stderr

# Set Java environment (required for AGeNT)
{java_env}

# Verify Java version
java -version

# Run AGeNT trim
{agent_path} trim \\
  -fq1 {fq1} \\
  -fq2 {fq2} \\
  -out {output_prefix} \\
  -adaptor {adaptor} \\
  -mbc {mbc}

echo "AGeNT trimming completed for {sample}"
"""

    job_script = job_script_template.format(
        sample=sample_name,
        memory=config['sge']['memory'],
        cores=config['sge']['cores'],
        wall_time=config['sge']['wall_time'],
        output_dir=output_dir,
        java_env=java_env,
        agent_path=config['tools']['agent_path'],
        fq1=fastq1,
        fq2=fastq2,
        output_prefix=output_prefix,
        adaptor=config['agent']['adaptor'],
        mbc=config['agent']['mbc_type']
    )

    job_script_path = os.path.join(script_dir, f"{sample_name}_agent_trim.sh")
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
        description='Generate and submit SGE jobs for AGeNT MBC trimming'
    )
    parser.add_argument(
        '--config',
        default='config/01_agent_trim_config.yaml',
        help='Path to configuration file (default: config/01_agent_trim_config.yaml)'
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

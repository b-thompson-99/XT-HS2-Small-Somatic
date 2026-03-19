#!/usr/bin/env python3
"""
Step 2a: CReaK Molecular Consensus Calling

Generates and submits SGE jobs to run AGeNT CReaK consensus calling at one or
more stringency levels. Groups reads by MBC and generates consensus sequences
for accurate low-frequency somatic variant detection.

Stringency levels:
    MS1    — HYBRID mode, min single strand consensus 1 (maximum sensitivity)
    MS2    — HYBRID mode, min single strand consensus 2 (balanced)
    DUPLEX — DUPLEX mode (maximum specificity, both strands in agreement required)

Input:
    - Directory containing sample subdirectories, each with a sorted BAM
      ending in *_sorted.bam (MBC tags must be present from steps 1a-1c)

Output:
    - Consensus BAMs in per-stringency subdirectories:
      {output_bam_dir}/MS1/{sample}_hybrid_MS1.bam
      {output_bam_dir}/MS2/{sample}_hybrid_MS2.bam
      {output_bam_dir}/DUPLEX/{sample}_duplex.bam

Requirements:
    - AGeNT (v3.1.2 or later)
    - Java 17+
    - SGE cluster with qsub
    - config_2a_agent_creak.yaml
"""

import os
import copy
import subprocess
import time
import yaml
import argparse
from pathlib import Path


def load_config(config_path):
    """Load configuration from YAML file."""
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)


def find_sorted_bam(sample_dir, bam_suffix):
    """
    Find sorted BAM file in a sample directory.

    Returns:
        Path to BAM file or None if not found.
    """
    bam_candidates = [f for f in os.listdir(sample_dir) if f.endswith(bam_suffix)]

    if len(bam_candidates) == 0:
        print(f"Warning: No BAM file found in {sample_dir}")
        return None
    elif len(bam_candidates) > 1:
        print(f"Warning: Multiple BAM files found in {sample_dir}")
        print(f"  Candidates: {bam_candidates}")
        print(f"  Using first match: {bam_candidates[0]}")

    return os.path.join(sample_dir, bam_candidates[0])


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


def get_output_suffix(consensus_mode, min_stringency):
    """Return output BAM suffix based on consensus parameters."""
    if consensus_mode.upper() == "DUPLEX":
        return "_duplex"
    elif consensus_mode.upper() == "HYBRID":
        return f"_hybrid_MS{min_stringency}"
    else:
        return "_consensus"


def create_qsub_script(sample_name, config, input_bam, output_dir, script_dir, stringency_label):
    """
    Create SGE job script for AGeNT CReaK consensus calling.

    Returns:
        Path to created job script, or None if output already exists.
    """
    output_suffix = get_output_suffix(
        config['creak']['consensus_mode'],
        config['creak']['min_stringency']
    )

    output_bam = os.path.join(output_dir, f"{sample_name}{output_suffix}.bam")

    if os.path.exists(output_bam) and os.path.getsize(output_bam) > 0:
        print(f"  ⚠ Output already exists: {os.path.basename(output_bam)} (skipping)")
        return None

    java_home = config['tools']['java_home']
    java_env = f"""export JAVA_HOME={java_home}
export PATH={java_home}/bin:$PATH"""

    creak_params = config['creak']
    creak_cmd = f"""{config['tools']['agent_path']} creak \\
  -MD {creak_params['max_distance']} \\
  -MS {creak_params['min_stringency']} \\
  -c={creak_params['consensus_mode']} \\
  -s={creak_params['sample_size']} \\
  -mm={creak_params['max_mismatches']} \\
  -mr={creak_params['max_read_length']} \\
  -mq={creak_params['min_mapping_quality']}"""

    if creak_params.get('memory_efficient', True):
        creak_cmd += " \\\n  --memory-efficient-mode"
    if creak_params.get('force_overwrite', True):
        creak_cmd += " \\\n  -F"
    if creak_params.get('filter_mode', True):
        creak_cmd += " \\\n  -f"
    if creak_params.get('retain_singletons', True):
        creak_cmd += " \\\n  -r"

    creak_cmd += f" \\\n  -o={output_bam} \\\n  {input_bam}"

    job_script_template = """#!/bin/bash
#$ -N creak_{sample}_{stringency}
#$ -l h_vmem={memory}
#$ -pe smp {cores}
#$ -l h_rt={wall_time}
#$ -cwd
#$ -o {log_dir}/{sample}_creak_{stringency}.stdout
#$ -e {log_dir}/{sample}_creak_{stringency}.stderr

# Set Java environment (required for AGeNT)
{java_env}

# Verify Java version
java -version

INPUT_BAM={input_bam}
OUTPUT_BAM={output_bam}

echo "Starting AGeNT CReaK for {sample} ({stringency})..."
echo "  Input BAM: $(basename $INPUT_BAM)"
echo "  Output BAM: $(basename $OUTPUT_BAM)"
echo "  Consensus mode: {consensus_mode}"
echo "  Minimum stringency: MS{min_stringency}"

# Run AGeNT CReaK
{creak_cmd}

# Verify output
if [ -f "$OUTPUT_BAM" ]; then
    SIZE=$(stat -f%z "$OUTPUT_BAM" 2>/dev/null || stat -c%s "$OUTPUT_BAM" 2>/dev/null)
    echo "✓ CReaK completed successfully"
    echo "  Output BAM: $(basename $OUTPUT_BAM)"
    echo "  File size: $SIZE bytes"
else
    echo "✗ ERROR: Output BAM not created"
    exit 1
fi
"""

    job_script = job_script_template.format(
        sample=sample_name,
        stringency=stringency_label,
        memory=config['sge']['memory'],
        cores=config['sge']['cores'],
        wall_time=config['sge']['wall_time'],
        log_dir=output_dir,
        java_env=java_env,
        input_bam=input_bam,
        output_bam=output_bam,
        consensus_mode=creak_params['consensus_mode'],
        min_stringency=creak_params['min_stringency'],
        creak_cmd=creak_cmd
    )

    job_script_path = os.path.join(script_dir, f"{sample_name}_creak_{stringency_label}.sh")
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


def run_stringency(stringency, config, dry_run):
    """
    Run CReaK pipeline for a single stringency level.

    Returns:
        List of job scripts generated.
    """
    cfg = copy.deepcopy(config)

    stringency_map = {
        'MS1':    {'mode': 'HYBRID',  'min_stringency': 1, 'suffix': 'MS1'},
        'MS2':    {'mode': 'HYBRID',  'min_stringency': 2, 'suffix': 'MS2'},
        'DUPLEX': {'mode': 'DUPLEX',  'min_stringency': 2, 'suffix': 'DUPLEX'}
    }

    stringency_params = stringency_map[stringency]
    cfg['creak']['consensus_mode'] = stringency_params['mode']
    cfg['creak']['min_stringency'] = stringency_params['min_stringency']

    base_output_dir = cfg['paths']['output_bam_dir']
    output_dir = os.path.join(base_output_dir, stringency_params['suffix'])
    script_dir = cfg['paths']['job_scripts_dir']

    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(script_dir, exist_ok=True)

    print(f"\n{'='*60}")
    print(f"Stringency:       {stringency}")
    print(f"Consensus Mode:   {cfg['creak']['consensus_mode']}")
    print(f"Min Stringency:   MS{cfg['creak']['min_stringency']}")
    print(f"Output Directory: {output_dir}")
    print(f"{'='*60}\n")

    input_dir = cfg['paths']['input_bam_dir']
    job_scripts = []

    print(f"Scanning for samples in: {input_dir}\n")

    for sample in sorted(os.listdir(input_dir)):
        sample_dir = os.path.join(input_dir, sample)

        if not os.path.isdir(sample_dir):
            continue

        print(f"Processing sample: {sample}")

        input_bam = find_sorted_bam(sample_dir, cfg['bam']['input_suffix'])

        if input_bam is None:
            print(f"  ⚠ Skipping {sample} - BAM file not found\n")
            continue

        print(f"  Input BAM: {os.path.basename(input_bam)}")

        job_script_path = create_qsub_script(
            sample_name=sample,
            config=cfg,
            input_bam=input_bam,
            output_dir=output_dir,
            script_dir=script_dir,
            stringency_label=stringency
        )

        if job_script_path is not None:
            job_scripts.append(job_script_path)
            print(f"  ✓ Job script created: {os.path.basename(job_script_path)}\n")
        else:
            print()

    return job_scripts


def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(
        description='Generate and submit AGeNT CReaK jobs for molecular consensus calling',
        epilog='Example: python 2a_agent_creak.py --config config_2a_agent_creak.yaml --stringency MS1 MS2 DUPLEX'
    )
    parser.add_argument(
        '--config',
        default='config/04_agent_creak_config.yaml',
        help='Path to configuration file (default: config/04_agent_creak_config.yaml)'
    )
    parser.add_argument(
        '--stringency',
        choices=['MS1', 'MS2', 'DUPLEX'],
        nargs='+',
        required=True,
        help='One or more stringency levels to run: MS1 MS2 DUPLEX'
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

    all_job_scripts = []

    for stringency in args.stringency:
        job_scripts = run_stringency(stringency, config, args.dry_run)
        all_job_scripts.extend(job_scripts)

    print(f"\n{'='*60}")
    print(f"Generated {len(all_job_scripts)} job scripts total")
    print(f"Stringencies processed: {', '.join(args.stringency)}")
    print(f"{'='*60}\n")

    if len(all_job_scripts) == 0:
        print("No jobs to submit. Exiting.")
        return

    if args.dry_run:
        print("DRY RUN: Jobs not submitted.")
        print("Job scripts saved to:", config['paths']['job_scripts_dir'])
    else:
        submit_jobs_in_batches(
            all_job_scripts,
            batch_size=config['submission']['batch_size'],
            pause_time=config['submission']['pause_between_batches']
        )
        print("\n✓ All batches submitted successfully")


if __name__ == "__main__":
    main()

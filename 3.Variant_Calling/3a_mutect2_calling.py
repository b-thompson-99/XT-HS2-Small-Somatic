#!/usr/bin/env python3
"""
Step 3a: Somatic Variant Calling with GATK Mutect2

Generates and submits SGE jobs to run Mutect2 in tumour-only (no matched sample) mode on consensus
BAMs from step 2b. Run separately for each stringency level (MS1, MS2, DUPLEX).

Input:
    - Directory containing consensus BAMs with read groups (With_rg/ from step 2b)

Output:
    - Compressed VCFs: {sample}.vcf.gz + .tbi index
    - Orientation bias files: {sample}_f1r2.tar.gz
    - Job logs: {sample}.log / .err

Requirements:
    - GATK (v4.6.1.0 or later)
    - Java 17+
    - Reference genome (hg38) with .fai index
    - Target BED file
    - SGE cluster with qsub
    - config_3a_mutect2.yaml
"""

import os
import subprocess
import argparse
import yaml
from pathlib import Path


def load_config(config_path):
    """Load configuration from YAML file."""
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)


def verify_dependencies(config):
    """
    Verify required tools and files exist.

    Raises:
        FileNotFoundError: If dependencies are missing.
    """
    gatk_jar = config['tools']['gatk_jar']
    reference = config['reference']['genome_fasta']
    bed_file = config['reference']['target_bed']
    java_home = config['tools']['java_home']

    if not os.path.exists(gatk_jar):
        raise FileNotFoundError(f"GATK JAR not found: {gatk_jar}")

    if not os.path.exists(reference):
        raise FileNotFoundError(f"Reference genome not found: {reference}")

    ref_index = reference + ".fai"
    if not os.path.exists(ref_index):
        raise FileNotFoundError(
            f"Reference index not found: {ref_index}\n"
            f"Run: samtools faidx {reference}"
        )

    if not os.path.exists(bed_file):
        raise FileNotFoundError(f"Target BED file not found: {bed_file}")

    if not os.path.exists(java_home):
        raise FileNotFoundError(f"Java installation not found: {java_home}")

    print(f"✓ GATK found: {gatk_jar}")
    print(f"✓ Reference genome found: {reference}")
    print(f"✓ Target BED found: {bed_file}")
    print(f"✓ Java found: {java_home}")


def find_bam_files(input_dir):
    """Find all BAM files in the input directory."""
    bam_files = []
    for filename in os.listdir(input_dir):
        if filename.endswith(".bam"):
            bam_files.append(os.path.join(input_dir, filename))
    return sorted(bam_files)


def create_mutect2_job(bam_file, config, output_dir, script_dir):
    """
    Create SGE job script for Mutect2 variant calling.

    Returns:
        Tuple of (job_script_path, sample_name), or (None, sample_name) if output exists.
    """
    sample_name = os.path.splitext(os.path.basename(bam_file))[0]

    vcf_file = os.path.join(output_dir, f"{sample_name}.vcf.gz")
    f1r2_tar = os.path.join(output_dir, f"{sample_name}_f1r2.tar.gz")
    log_file = os.path.join(output_dir, f"{sample_name}.log")
    err_file = os.path.join(output_dir, f"{sample_name}.err")

    if os.path.exists(vcf_file) and os.path.exists(vcf_file + ".tbi"):
        if os.path.getsize(vcf_file) > 0:
            print(f"  ⚠ Output already exists: {sample_name}.vcf.gz (skipping)")
            return None, sample_name

    job_script_template = """#!/bin/bash
#$ -N mutect2_{sample}
#$ -l h_vmem={memory}
#$ -pe smp {cores}
#$ -l h_rt={wall_time}
#$ -cwd
#$ -o {log_file}
#$ -e {err_file}
#$ -V

# Set Java environment
export JAVA_HOME={java_home}
export PATH=$JAVA_HOME/bin:$PATH

# Verify Java version
java -version

INPUT_BAM={bam_file}
OUTPUT_VCF={vcf_file}
F1R2_TAR={f1r2_tar}
REFERENCE={reference}
BED_FILE={bed_file}

echo "Starting Mutect2 for {sample}..."
echo "  Input BAM: $(basename $INPUT_BAM)"
echo "  Output VCF: $(basename $OUTPUT_VCF)"
echo "  Target regions: $(basename $BED_FILE)"

# Run GATK Mutect2
java -Xmx{java_memory} -jar {gatk_jar} Mutect2 \\
    -R $REFERENCE \\
    -I $INPUT_BAM \\
    -O $OUTPUT_VCF \\
    -L $BED_FILE \\
    --tumor-sample {sample} \\
    --f1r2-tar-gz $F1R2_TAR

# Index the output VCF
echo "Indexing VCF..."
java -jar {gatk_jar} IndexFeatureFile -I $OUTPUT_VCF

# Verify output
if [ -f "$OUTPUT_VCF" ] && [ -f "$OUTPUT_VCF.tbi" ]; then
    SIZE=$(stat -f%z "$OUTPUT_VCF" 2>/dev/null || stat -c%s "$OUTPUT_VCF" 2>/dev/null)
    VARIANTS=$(zcat "$OUTPUT_VCF" | grep -v "^#" | wc -l)
    echo "✓ Mutect2 completed successfully"
    echo "  Output VCF: $(basename $OUTPUT_VCF)"
    echo "  File size: $SIZE bytes"
    echo "  Variant count: $VARIANTS"
else
    echo "✗ ERROR: Output VCF or index not created"
    exit 1
fi
"""

    total_memory_gb = int(config['sge']['memory'].replace('G', ''))
    java_memory = f"{total_memory_gb - 1}G"

    job_script = job_script_template.format(
        sample=sample_name,
        memory=config['sge']['memory'],
        cores=config['sge']['cores'],
        wall_time=config['sge']['wall_time'],
        log_file=log_file,
        err_file=err_file,
        java_home=config['tools']['java_home'],
        bam_file=bam_file,
        vcf_file=vcf_file,
        f1r2_tar=f1r2_tar,
        reference=config['reference']['genome_fasta'],
        bed_file=config['reference']['target_bed'],
        gatk_jar=config['tools']['gatk_jar'],
        java_memory=java_memory
    )

    job_script_path = os.path.join(script_dir, f"{sample_name}_mutect2.sh")
    with open(job_script_path, "w") as f:
        f.write(job_script)

    return job_script_path, sample_name


def create_completion_check_job(output_dir, script_dir, mutect_job_ids):
    """
    Create a job to verify all Mutect2 outputs are complete.
    Runs after all Mutect2 jobs finish via SGE job dependencies.

    Returns:
        Path to completion check job script.
    """
    check_script = os.path.join(script_dir, "check_mutect2_completion.sh")
    check_log = os.path.join(output_dir, "check_mutect2_completion.log")
    check_err = os.path.join(output_dir, "check_mutect2_completion.err")

    check_script_template = """#!/bin/bash
#$ -N check_mutect2
#$ -l h_rt=00:15:00
#$ -cwd
#$ -o {check_log}
#$ -e {check_err}
#$ -V

echo "=========================================="
echo "Checking Mutect2 Completion"
echo "=========================================="
echo "Output directory: {output_dir}"
echo ""

incomplete=0
total=0

for vcf in {output_dir}/*.vcf.gz; do
    [ -e "$vcf" ] || continue
    total=$((total + 1))
    sample=$(basename "$vcf" .vcf.gz)

    if [ ! -s "$vcf" ]; then
        echo "❌ Missing or empty: $sample.vcf.gz"
        incomplete=$((incomplete + 1))
        continue
    fi

    if [ ! -f "$vcf.tbi" ]; then
        echo "❌ Index missing: $sample.vcf.gz.tbi"
        incomplete=$((incomplete + 1))
        continue
    fi

    variant_count=$(zcat "$vcf" | grep -v "^#" | wc -l)
    echo "✅ $sample: $variant_count variants"
done

echo ""
echo "=========================================="
echo "Summary"
echo "=========================================="
echo "Total samples: $total"
echo "Complete: $((total - incomplete))"
echo "Incomplete: $incomplete"

if [ "$incomplete" -eq 0 ]; then
    echo ""
    echo "✅✅ All Mutect2 outputs complete!"
    exit 0
else
    echo ""
    echo "⚠️ Some outputs are incomplete. Check logs above."
    exit 1
fi
"""

    check_script_content = check_script_template.format(
        output_dir=output_dir,
        check_log=check_log,
        check_err=check_err
    )

    with open(check_script, "w") as f:
        f.write(check_script_content)

    return check_script


def submit_jobs(job_scripts, output_dir, script_dir, dry_run=False):
    """
    Submit Mutect2 jobs and a completion check job.

    Returns:
        List of submitted job IDs.
    """
    mutect_job_ids = []

    print(f"\n{'='*60}")
    print(f"Submitting {len(job_scripts)} Mutect2 jobs")
    print(f"{'='*60}\n")

    for job_script in job_scripts:
        if dry_run:
            print(f"  [DRY RUN] Would submit: {os.path.basename(job_script)}")
            mutect_job_ids.append("12345")
        else:
            result = subprocess.run(
                ["qsub", job_script],
                capture_output=True,
                text=True
            )

            if result.returncode == 0:
                job_id = result.stdout.strip().split()[2]
                mutect_job_ids.append(job_id)
                print(f"  ✓ Submitted: {os.path.basename(job_script)}")
                print(f"    Job ID: {job_id}")
            else:
                print(f"  ✗ Failed: {os.path.basename(job_script)}")
                print(f"    Error: {result.stderr.strip()}")

    if mutect_job_ids:
        check_script = create_completion_check_job(output_dir, script_dir, mutect_job_ids)

        if dry_run:
            print(f"\n  [DRY RUN] Would submit completion check job")
            print(f"    Dependencies: {','.join(mutect_job_ids)}")
        else:
            hold_str = ",".join(mutect_job_ids)
            result = subprocess.run(
                ["qsub", "-hold_jid", hold_str, check_script],
                capture_output=True,
                text=True
            )

            if result.returncode == 0:
                print(f"\n✓ Submitted completion check job")
                print(f"  Dependencies: {hold_str}")
            else:
                print(f"\n⚠️ Failed to submit completion check job")

    return mutect_job_ids


def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(
        description='Generate and submit GATK Mutect2 jobs for somatic variant calling',
        epilog='Example: python 3a_mutect2_calling.py --config config_3a_mutect2.yaml --input /path/to/creak_output/MS1/With_rg'
    )
    parser.add_argument(
        '--config',
        default='config/05_mutect2_config.yaml',
        help='Path to configuration file (default: config/05_mutect2_config.yaml)'
    )
    parser.add_argument(
        '--input',
        help='Input BAM directory — overrides config. Use With_rg/ subdirectory from step 2b.'
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Generate job scripts but do not submit them'
    )

    args = parser.parse_args()

    config = load_config(args.config)

    if args.input:
        config['paths']['input_bam_dir'] = os.path.abspath(args.input)
        config['paths']['output_vcf_dir'] = os.path.join(args.input, "Mutect_output")

    print("Checking dependencies...")
    verify_dependencies(config)
    print()

    input_dir = config['paths']['input_bam_dir']
    output_dir = config['paths']['output_vcf_dir']
    script_dir = config['paths']['job_scripts_dir']

    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(script_dir, exist_ok=True)

    bam_files = find_bam_files(input_dir)

    if len(bam_files) == 0:
        print(f"ERROR: No BAM files found in {input_dir}")
        return

    print(f"Found {len(bam_files)} BAM files in: {input_dir}\n")

    job_scripts = []

    for bam_file in bam_files:
        sample_name = os.path.splitext(os.path.basename(bam_file))[0]
        print(f"Processing: {sample_name}")

        job_script_path, _ = create_mutect2_job(
            bam_file=bam_file,
            config=config,
            output_dir=output_dir,
            script_dir=script_dir
        )

        if job_script_path is not None:
            job_scripts.append(job_script_path)
            print(f"  ✓ Job script created: {os.path.basename(job_script_path)}\n")
        else:
            print()

    print(f"{'='*60}")
    print(f"Generated {len(job_scripts)} job scripts")
    print(f"{'='*60}\n")

    if len(job_scripts) == 0:
        print("No jobs to submit. All samples already processed.")
        return

    if args.dry_run:
        print("DRY RUN: Jobs not submitted.")
        print(f"Job scripts saved to: {script_dir}")
    else:
        job_ids = submit_jobs(job_scripts, output_dir, script_dir, dry_run=False)
        print(f"\n✓ Submitted {len(job_ids)} Mutect2 jobs")


if __name__ == "__main__":
    main()

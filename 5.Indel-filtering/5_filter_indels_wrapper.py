#!/usr/bin/env python3
"""
Step 5: Indel Filtering Pipeline Wrapper

Runs the four-stage indel filtering pipeline sequentially:
1. Allele fraction (AF) filtering
2. Strand bias (SOR) filtering
3. Read quality filtering
4. RepeatMasker exclusion

Input:
    - Duplex indel VCFs from step 3b (Candidates/Indels/)

Usage:
    python 5_filter_indels_wrapper.py \
        --input /path/to/Candidates/Indels \
        --repeat-bed /path/to/hg38_repeats.bed
"""

import os
import sys
import argparse
import subprocess


def run_filter_step(script_name, input_dir, args_dict, step_name, dry_run=False):
    """
    Run a single filtering step.

    Returns:
        True if successful, False otherwise.
    """
    cmd = ['python3', script_name, '--input', input_dir]
    for arg_name, arg_value in args_dict.items():
        cmd.extend([f'--{arg_name}', str(arg_value)])

    if dry_run:
        print(f"  [DRY RUN] Would execute: {' '.join(cmd)}\n")
        return True

    print(f"Running {step_name}...")
    print(f"  Command: {' '.join(cmd)}")
    print()

    result = subprocess.run(cmd, capture_output=False, text=True)

    if result.returncode != 0:
        print(f"\nERROR: {step_name} failed with exit code {result.returncode}")
        return False

    print(f"\n✓ {step_name} complete\n")
    return True


def find_vcf_files(directory):
    """Count VCF files in directory."""
    if not os.path.isdir(directory):
        return 0
    return sum(1 for f in os.listdir(directory) if f.endswith('.vcf') or f.endswith('.vcf.gz'))


def main():
    parser = argparse.ArgumentParser(
        description='Run complete 4-stage indel filtering pipeline',
        epilog='Example: python 5_filter_indels_wrapper.py --input /path/to/Candidates/Indels --repeat-bed /path/to/hg38_repeats.bed'
    )
    parser.add_argument('--input', required=True,
                        help='Input directory containing unfiltered indel VCF files from step 3b')
    parser.add_argument('--repeat-bed', required=True,
                        help='Path to RepeatMasker BED file (e.g. hg38_repeats.bed)')
    parser.add_argument('--af-threshold', type=float, default=0.3,
                        help='AF threshold for step 5a (default: 0.3)')
    parser.add_argument('--sor-threshold', type=float, default=3.0,
                        help='SOR threshold for step 5b (default: 3.0)')
    parser.add_argument('--min-mmq', type=int, default=40,
                        help='Minimum mapping quality for step 5c (default: 40)')
    parser.add_argument('--max-ecnt', type=int, default=2,
                        help='Maximum event count for step 5c (default: 2)')
    parser.add_argument('--min-mpos', type=int, default=5,
                        help='Minimum distance from read end for step 5c (default: 5)')
    parser.add_argument('--dry-run', action='store_true',
                        help='Show what would be executed without running')

    args = parser.parse_args()

    if not os.path.isdir(args.input):
        print(f"ERROR: Input directory not found: {args.input}")
        return 1

    if not os.path.exists(args.repeat_bed):
        print(f"ERROR: RepeatMasker BED file not found: {args.repeat_bed}")
        return 1

    vcf_count = find_vcf_files(args.input)
    if vcf_count == 0:
        print(f"ERROR: No VCF files found in {args.input}")
        return 1

    print(f"\n{'='*60}")
    print(f"Indel Filtering Pipeline (4 Stages)")
    print(f"{'='*60}\n")
    print(f"Input directory: {args.input}")
    print(f"Found {vcf_count} VCF files")
    print(f"RepeatMasker BED: {args.repeat_bed}")

    if args.dry_run:
        print("\nDRY RUN MODE: No files will be modified\n")

    print(f"\nFilter parameters:")
    print(f"  AF threshold:  ≤{args.af_threshold}")
    print(f"  SOR threshold: ≤{args.sor_threshold}")
    print(f"  MMQ threshold: ≥{args.min_mmq}")
    print(f"  ECNT threshold: ≤{args.max_ecnt}")
    print(f"  MPOS threshold: ≥{args.min_mpos}")
    print()

    script_dir = os.path.dirname(os.path.abspath(__file__))

    input_5a = args.input
    input_5b = os.path.join(args.input, "AF_filtered")
    input_5c = os.path.join(args.input, "AF_filtered", "SB_filtered")
    input_5d = os.path.join(args.input, "AF_filtered", "SB_filtered", "Read_filtered")
    output_final = os.path.join(input_5d, "Repeats_removed")

    print(f"{'='*60}")
    print(f"Step 5a: Allele Fraction Filtering")
    print(f"{'='*60}\n")

    if not run_filter_step(
        script_name=os.path.join(script_dir, '5a_indel_filter_AF.py'),
        input_dir=input_5a,
        args_dict={'af-threshold': args.af_threshold},
        step_name='Step 5a (AF filtering)',
        dry_run=args.dry_run
    ):
        return 1

    print(f"{'='*60}")
    print(f"Step 5b: Strand Bias (SOR) Filtering")
    print(f"{'='*60}\n")

    if not run_filter_step(
        script_name=os.path.join(script_dir, '5b_indel_filter_SOR.py'),
        input_dir=input_5b,
        args_dict={'sor-threshold': args.sor_threshold},
        step_name='Step 5b (Strand bias filtering)',
        dry_run=args.dry_run
    ):
        return 1

    print(f"{'='*60}")
    print(f"Step 5c: Read Quality Filtering")
    print(f"{'='*60}\n")

    if not run_filter_step(
        script_name=os.path.join(script_dir, '5c_indel_filter_read_quality.py'),
        input_dir=input_5c,
        args_dict={
            'min-mmq': args.min_mmq,
            'max-ecnt': args.max_ecnt,
            'min-mpos': args.min_mpos
        },
        step_name='Step 5c (Read quality filtering)',
        dry_run=args.dry_run
    ):
        return 1

    print(f"{'='*60}")
    print(f"Step 5d: RepeatMasker Exclusion")
    print(f"{'='*60}\n")

    if not run_filter_step(
        script_name=os.path.join(script_dir, '5d_indel_exclude_repeats.py'),
        input_dir=input_5d,
        args_dict={'repeat-bed': args.repeat_bed},
        step_name='Step 5d (RepeatMasker exclusion)',
        dry_run=args.dry_run
    ):
        return 1

    print(f"\n{'='*60}")
    print(f"Indel Filtering Pipeline Complete!")
    print(f"{'='*60}\n")

    if not args.dry_run:
        final_vcf_count = find_vcf_files(output_final)
        print(f"✓ All 4 filtering stages completed successfully")
        print(f"  Input VCFs:  {vcf_count}")
        print(f"  Output VCFs: {final_vcf_count}")
        print(f"\nFinal filtered VCFs saved in:")
        print(f"  {output_final}/")
    else:
        print("DRY RUN: No files were modified")

    print(f"{'='*60}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())

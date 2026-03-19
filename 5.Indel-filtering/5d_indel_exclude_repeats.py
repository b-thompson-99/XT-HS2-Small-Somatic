#!/usr/bin/env python3
"""
Step 5d: Indel RepeatMasker Exclusion

Removes indels overlapping RepeatMasker regions using bedtools intersect.
Indels in repetitive regions (homopolymers, simple repeats, SINEs, LINEs etc.)
are highly error-prone even with Duplex consensus.

Input:
    - Read-filtered indel VCFs from step 5c (Read_filtered/)
    - RepeatMasker BED file (hg38_repeats.bed)

Output:
    - Filtered VCFs in: {input_dir}/Repeats_removed/
    - filtering_stats.txt

Requirements:
    - bedtools

Usage:
    python 5d_indel_exclude_repeats.py \
        --input /path/to/Candidates/Indels/AF_filtered/SB_filtered/Read_filtered \
        --repeat-bed /path/to/hg38_repeats.bed
"""

import os
import sys
import argparse
import subprocess
import gzip
from pathlib import Path


def run_command(cmd):
    """Run shell command and handle errors."""
    try:
        result = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
        return result
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {cmd}")
        print(f"Error: {e.stderr}")
        sys.exit(1)


def count_variants_in_vcf(vcf_file):
    """Count non-header lines in VCF file."""
    count = 0
    opener = gzip.open(vcf_file, 'rt') if vcf_file.endswith('.gz') else open(vcf_file, 'r')
    with opener as f:
        for line in f:
            if not line.startswith('#'):
                count += 1
    return count


def filter_vcf(input_vcf, output_vcf, repeat_bed):
    """
    Filter VCF to remove variants overlapping repeat regions using bedtools.

    Returns:
        Tuple of (variants_before, variants_after).
    """
    variants_before = count_variants_in_vcf(input_vcf)

    if input_vcf.endswith('.gz'):
        temp_vcf = f"temp_{os.path.basename(input_vcf).replace('.gz', '')}"
        run_command(f"gunzip -c {input_vcf} > {temp_vcf}")
        actual_input = temp_vcf
    else:
        actual_input = input_vcf

    cmd = f"bedtools intersect -header -v -a {actual_input} -b {repeat_bed} > {output_vcf}"
    run_command(cmd)

    if input_vcf.endswith('.gz'):
        run_command(f"rm {temp_vcf}")

    variants_after = count_variants_in_vcf(output_vcf)
    return variants_before, variants_after


def find_vcf_files(input_dir):
    """Find all VCF files in directory."""
    vcf_files = []
    for filename in os.listdir(input_dir):
        if filename.endswith(".vcf") or filename.endswith(".vcf.gz"):
            vcf_files.append(os.path.join(input_dir, filename))
    return sorted(vcf_files)


def process_directory(input_dir, repeat_bed, dry_run=False):
    """
    Process all VCF files in a directory.

    Returns:
        Dictionary of statistics per sample.
    """
    output_dir = os.path.join(input_dir, "Repeats_removed")

    if not dry_run:
        os.makedirs(output_dir, exist_ok=True)

    vcf_files = find_vcf_files(input_dir)

    if len(vcf_files) == 0:
        print(f"ERROR: No VCF files found in {input_dir}")
        return {}

    print(f"Found {len(vcf_files)} VCF files in: {input_dir}")
    print(f"RepeatMasker BED: {repeat_bed}")
    print(f"Output directory: {output_dir}\n")

    stats = {}

    for vcf_file in vcf_files:
        sample_name = os.path.basename(vcf_file).replace('.vcf.gz', '').replace('.vcf', '')
        output_name = os.path.basename(vcf_file).replace('.vcf.gz', '.vcf')
        output_vcf = os.path.join(output_dir, output_name)

        print(f"Processing: {sample_name}")

        if dry_run:
            print(f"  [DRY RUN] Would filter: {os.path.basename(vcf_file)}")
            print(f"  [DRY RUN] Output: {output_name}\n")
            continue

        before, after = filter_vcf(vcf_file, output_vcf, repeat_bed)
        removed = before - after

        stats[sample_name] = {'before': before, 'after': after, 'removed': removed}

        print(f"  Before: {before} variants")
        print(f"  After:  {after} variants")
        print(f"  Removed: {removed} variants ({removed/before*100:.1f}%)")
        print()

    return stats


def print_summary(stats):
    """Print summary statistics."""
    if not stats:
        return
    total_before = sum(s['before'] for s in stats.values())
    total_after = sum(s['after'] for s in stats.values())
    total_removed = sum(s['removed'] for s in stats.values())

    print(f"{'='*60}")
    print(f"SUMMARY: Indel RepeatMasker Exclusion")
    print(f"{'='*60}")
    print(f"Samples processed: {len(stats)}")
    print(f"Total variants before: {total_before}")
    print(f"Total variants after:  {total_after}")
    print(f"Total variants removed: {total_removed} ({total_removed/total_before*100:.1f}%)")
    print(f"{'='*60}\n")


def write_stats_file(stats, output_dir):
    """Write filtering statistics to text file."""
    stats_file = os.path.join(output_dir, "filtering_stats.txt")
    total_before = sum(s['before'] for s in stats.values())
    total_after = sum(s['after'] for s in stats.values())
    total_removed = sum(s['removed'] for s in stats.values())

    with open(stats_file, 'w') as f:
        f.write("="*60 + "\n")
        f.write("Indel RepeatMasker Exclusion Statistics\n")
        f.write("="*60 + "\n\n")
        for sample_name, s in sorted(stats.items()):
            f.write(f"{sample_name}:\n")
            f.write(f"  Input variants:  {s['before']}\n")
            f.write(f"  Output variants: {s['after']}\n")
            f.write(f"  Removed:         {s['removed']} ({s['removed']/s['before']*100:.1f}%)\n\n")
        f.write("="*60 + "\n")
        f.write("TOTAL:\n")
        f.write(f"  Samples processed: {len(stats)}\n")
        f.write(f"  Input variants:    {total_before}\n")
        f.write(f"  Output variants:   {total_after}\n")
        f.write(f"  Removed:           {total_removed} ({total_removed/total_before*100:.1f}%)\n")
        f.write("="*60 + "\n")

    print(f"✓ Statistics saved to: {stats_file}")


def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(
        description='Filter indel VCF files to remove variants in RepeatMasker regions',
        epilog='Example: python 5d_indel_exclude_repeats.py --input /path/to/Read_filtered --repeat-bed /path/to/hg38_repeats.bed'
    )
    parser.add_argument('--input', required=True,
                        help='Input directory containing read-filtered indel VCF files')
    parser.add_argument('--repeat-bed', required=True,
                        help='Path to RepeatMasker BED file (e.g. hg38_repeats.bed)')
    parser.add_argument('--dry-run', action='store_true',
                        help='Show what would be processed without filtering')

    args = parser.parse_args()

    if not os.path.isdir(args.input):
        print(f"ERROR: Input directory not found: {args.input}")
        return

    if not os.path.exists(args.repeat_bed):
        print(f"ERROR: RepeatMasker BED file not found: {args.repeat_bed}")
        print(f"\nDownload from UCSC Table Browser:")
        print(f"  Assembly: hg38, Group: Repeats, Track: RepeatMasker, Output format: BED")
        return

    result = subprocess.run(['bedtools', '--version'], capture_output=True)
    if result.returncode != 0:
        print("ERROR: bedtools not found — install with: conda install -c bioconda bedtools")
        return

    print(f"\n{'='*60}")
    print(f"Step 5d: Indel RepeatMasker Exclusion")
    print(f"{'='*60}\n")

    stats = process_directory(args.input, args.repeat_bed, args.dry_run)

    if not args.dry_run:
        print_summary(stats)
        output_dir = os.path.join(args.input, "Repeats_removed")
        write_stats_file(stats, output_dir)
        print("✓ Indel RepeatMasker exclusion complete!")
        print(f"  Filtered VCFs saved in: {output_dir}/\n")
        print("These are your final candidate indels (pre-IGV validation)")
    else:
        print("DRY RUN: No files were modified.\n")


if __name__ == "__main__":
    main()

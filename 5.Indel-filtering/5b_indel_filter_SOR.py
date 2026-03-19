#!/usr/bin/env python3
"""
Step 5b: Indel Strand Bias (SOR) Filtering

Removes strand-biased indels using the Strand Odds Ratio (SOR) metric computed
from the SB field in Mutect2 VCFs. SOR > 3 indicates likely artefact.

Input:
    - AF-filtered indel VCFs from step 5a (AF_filtered/)

Output:
    - Filtered VCFs in: {input_dir}/SB_filtered/
    - filtering_stats.txt

Usage:
    python 5b_indel_filter_SOR.py --input /path/to/Candidates/Indels/AF_filtered
    python 5b_indel_filter_SOR.py --input /path/to/Candidates/Indels/AF_filtered --sor-threshold 3
"""

import os
import argparse
import gzip
import math


def parse_sb_from_format(line):
    """
    Parse SB field from VCF FORMAT.

    Returns:
        List of 4 values [ref_fwd, ref_rev, alt_fwd, alt_rev] or None.
    """
    columns = line.strip().split('\t')
    if len(columns) < 10:
        return None

    format_keys = columns[8].split(':')
    sample_values = columns[9].split(':')

    if len(format_keys) != len(sample_values):
        return None

    try:
        sb_index = format_keys.index('SB')
        sb_parts = sample_values[sb_index].split(',')
        if len(sb_parts) == 4:
            return [int(x) for x in sb_parts]
    except (ValueError, IndexError):
        pass

    return None


def compute_sor(sb_values):
    """
    Compute Strand Odds Ratio (SOR) from SB field values.

    Returns:
        SOR value as float.
    """
    ref_fwd, ref_rev, alt_fwd, alt_rev = sb_values

    ref_fwd += 1; ref_rev += 1; alt_fwd += 1; alt_rev += 1

    ref_ratio = ref_fwd / ref_rev
    alt_ratio = alt_fwd / alt_rev

    sor = ref_ratio / alt_ratio if ref_ratio > alt_ratio else alt_ratio / ref_ratio
    return math.log(sor)


def filter_vcf_by_strand_bias(input_vcf, output_vcf, sor_threshold):
    """
    Filter VCF by strand odds ratio threshold.

    Returns:
        Tuple of (variants_before, variants_after, variants_removed, variants_no_sb).
    """
    variants_before = 0
    variants_after = 0
    variants_no_sb = 0

    opener = gzip.open(input_vcf, 'rt') if input_vcf.endswith('.gz') else open(input_vcf, 'r')

    with opener as infile, open(output_vcf, 'w') as outfile:
        for line in infile:
            if line.startswith("##") or line.startswith("#"):
                outfile.write(line)
                continue

            variants_before += 1
            sb_values = parse_sb_from_format(line)

            if sb_values is None:
                outfile.write(line)
                variants_after += 1
                variants_no_sb += 1
                continue

            sor = compute_sor(sb_values)

            if sor <= sor_threshold:
                outfile.write(line)
                variants_after += 1

    variants_removed = variants_before - variants_after
    return variants_before, variants_after, variants_removed, variants_no_sb


def find_vcf_files(input_dir):
    """Find all VCF files in directory."""
    vcf_files = []
    for filename in os.listdir(input_dir):
        if filename.endswith(".vcf") or filename.endswith(".vcf.gz"):
            vcf_files.append(os.path.join(input_dir, filename))
    return sorted(vcf_files)


def process_directory(input_dir, sor_threshold, dry_run=False):
    """
    Process all VCF files in a directory.

    Returns:
        Dictionary of statistics per sample.
    """
    output_dir = os.path.join(input_dir, "SB_filtered")

    if not dry_run:
        os.makedirs(output_dir, exist_ok=True)

    vcf_files = find_vcf_files(input_dir)

    if len(vcf_files) == 0:
        print(f"ERROR: No VCF files found in {input_dir}")
        return {}

    print(f"Found {len(vcf_files)} VCF files in: {input_dir}")
    print(f"SOR threshold: >{sor_threshold} will be removed (keeping ≤{sor_threshold})")
    print(f"Output directory: {output_dir}\n")

    stats = {}

    for vcf_file in vcf_files:
        sample_name = os.path.basename(vcf_file).replace('.vcf.gz', '').replace('.vcf', '')
        output_vcf = os.path.join(output_dir, os.path.basename(vcf_file).replace('.vcf.gz', '.vcf'))

        print(f"Processing: {sample_name}")

        if dry_run:
            print(f"  [DRY RUN] Would filter: {os.path.basename(vcf_file)}")
            print(f"  [DRY RUN] Output: {os.path.basename(output_vcf)}\n")
            continue

        before, after, removed, no_sb = filter_vcf_by_strand_bias(vcf_file, output_vcf, sor_threshold)

        stats[sample_name] = {'before': before, 'after': after, 'removed': removed, 'no_sb': no_sb}

        print(f"  Before: {before} variants")
        print(f"  After:  {after} variants")
        print(f"  Removed: {removed} variants ({removed/before*100:.1f}%)")
        if no_sb > 0:
            print(f"  No SB field: {no_sb} variants (kept)")
        print()

    return stats


def print_summary(stats):
    """Print summary statistics."""
    if not stats:
        return
    total_before = sum(s['before'] for s in stats.values())
    total_after = sum(s['after'] for s in stats.values())
    total_removed = sum(s['removed'] for s in stats.values())
    total_no_sb = sum(s['no_sb'] for s in stats.values())

    print(f"{'='*60}")
    print(f"SUMMARY: Indel Strand Bias (SOR) Filtering")
    print(f"{'='*60}")
    print(f"Samples processed: {len(stats)}")
    print(f"Total variants before: {total_before}")
    print(f"Total variants after:  {total_after}")
    print(f"Total variants removed: {total_removed} ({total_removed/total_before*100:.1f}%)")
    if total_no_sb > 0:
        print(f"Variants without SB field: {total_no_sb} (kept)")
    print(f"{'='*60}\n")


def write_stats_file(stats, output_dir):
    """Write filtering statistics to text file."""
    stats_file = os.path.join(output_dir, "filtering_stats.txt")
    total_before = sum(s['before'] for s in stats.values())
    total_after = sum(s['after'] for s in stats.values())
    total_removed = sum(s['removed'] for s in stats.values())
    total_no_sb = sum(s['no_sb'] for s in stats.values())

    with open(stats_file, 'w') as f:
        f.write("="*60 + "\n")
        f.write("Indel Strand Bias (SOR) Filtering Statistics\n")
        f.write("="*60 + "\n\n")
        for sample_name, s in sorted(stats.items()):
            f.write(f"{sample_name}:\n")
            f.write(f"  Input variants:  {s['before']}\n")
            f.write(f"  Output variants: {s['after']}\n")
            f.write(f"  Removed:         {s['removed']} ({s['removed']/s['before']*100:.1f}%)\n")
            if s['no_sb'] > 0:
                f.write(f"  No SB field:     {s['no_sb']}\n")
            f.write("\n")
        f.write("="*60 + "\n")
        f.write("TOTAL:\n")
        f.write(f"  Samples processed: {len(stats)}\n")
        f.write(f"  Input variants:    {total_before}\n")
        f.write(f"  Output variants:   {total_after}\n")
        f.write(f"  Removed:           {total_removed} ({total_removed/total_before*100:.1f}%)\n")
        if total_no_sb > 0:
            f.write(f"  No SB field:       {total_no_sb}\n")
        f.write("="*60 + "\n")

    print(f"✓ Statistics saved to: {stats_file}")


def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(
        description='Filter indel VCF files by strand bias (SOR)',
        epilog='Example: python 5b_indel_filter_SOR.py --input /path/to/Candidates/Indels/AF_filtered'
    )
    parser.add_argument('--input', required=True,
                        help='Input directory containing AF-filtered indel VCF files')
    parser.add_argument('--sor-threshold', type=float, default=3.0,
                        help='Maximum SOR to keep (default: 3.0)')
    parser.add_argument('--dry-run', action='store_true',
                        help='Show what would be processed without filtering')

    args = parser.parse_args()

    if not os.path.isdir(args.input):
        print(f"ERROR: Input directory not found: {args.input}")
        return

    if args.sor_threshold < 0:
        print(f"ERROR: SOR threshold must be positive (got {args.sor_threshold})")
        return

    print(f"\n{'='*60}")
    print(f"Step 5b: Indel Strand Bias (SOR) Filtering")
    print(f"{'='*60}\n")

    stats = process_directory(args.input, args.sor_threshold, args.dry_run)

    if not args.dry_run:
        print_summary(stats)
        output_dir = os.path.join(args.input, "SB_filtered")
        write_stats_file(stats, output_dir)
        print("✓ Indel strand bias filtering complete!")
        print(f"  Filtered VCFs saved in: {output_dir}/\n")
    else:
        print("DRY RUN: No files were modified.\n")


if __name__ == "__main__":
    main()

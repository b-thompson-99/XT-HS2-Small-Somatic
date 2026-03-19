#!/usr/bin/env python3
"""
Step 4a: Allele Fraction (AF) Filtering

Removes high allele fraction variants (AF > threshold) likely to be germline.
Applied to VCFs from all stringency levels (MS1, MS2, Duplex).

Input:
    - Directory containing SNV VCF files from step 3b

Output:
    - Filtered VCFs in: {input_dir}/AF_filtered/
    - filtering_stats.txt

Usage:
    python 4a_filter_AF.py --input /path/to/Candidates/SNVs-MS1 --af-threshold 0.3
"""

import os
import argparse
from pathlib import Path


def filter_vcf_by_af(input_vcf, output_vcf, af_threshold):
    """
    Filter VCF by allele fraction threshold.

    Returns:
        Tuple of (variants_before, variants_after, variants_removed).
    """
    variants_before = 0
    variants_after = 0

    with open(input_vcf, 'r') as infile, open(output_vcf, 'w') as outfile:
        for line in infile:
            if line.startswith("##") or line.startswith("#"):
                outfile.write(line)
                continue

            variants_before += 1
            columns = line.strip().split("\t")
            format_field = columns[8]
            sample_data = columns[9:]

            if 'AF' not in format_field.split(":"):
                outfile.write(line)
                variants_after += 1
                continue

            af_index = format_field.split(":").index('AF')

            keep_variant = True
            for sample in sample_data:
                af_value = sample.split(":")[af_index]
                try:
                    af_value = float(af_value)
                    if af_value > af_threshold:
                        keep_variant = False
                        break
                except (ValueError, IndexError):
                    continue

            if keep_variant:
                outfile.write(line)
                variants_after += 1

    variants_removed = variants_before - variants_after
    return variants_before, variants_after, variants_removed


def find_vcf_files(input_dir):
    """Find all VCF files in directory."""
    vcf_files = []
    for filename in os.listdir(input_dir):
        if filename.endswith(".vcf") or filename.endswith(".vcf.gz"):
            vcf_files.append(os.path.join(input_dir, filename))
    return sorted(vcf_files)


def process_directory(input_dir, af_threshold, dry_run=False):
    """
    Process all VCF files in a directory.

    Returns:
        Dictionary of statistics per sample.
    """
    output_dir = os.path.join(input_dir, "AF_filtered")

    if not dry_run:
        os.makedirs(output_dir, exist_ok=True)

    vcf_files = find_vcf_files(input_dir)

    if len(vcf_files) == 0:
        print(f"ERROR: No VCF files found in {input_dir}")
        return {}

    print(f"Found {len(vcf_files)} VCF files in: {input_dir}")
    print(f"AF threshold: >{af_threshold} will be removed (keeping ≤{af_threshold})")
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

        if vcf_file.endswith('.gz'):
            import gzip, shutil
            temp_vcf = vcf_file.replace('.gz', '.tmp')
            with gzip.open(vcf_file, 'rb') as f_in, open(temp_vcf, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
            vcf_file = temp_vcf

        before, after, removed = filter_vcf_by_af(vcf_file, output_vcf, af_threshold)

        if vcf_file.endswith('.tmp'):
            os.remove(vcf_file)

        stats[sample_name] = {
            'before': before,
            'after': after,
            'removed': removed,
            'percent_kept': (after / before * 100) if before > 0 else 0
        }

        print(f"  Before: {before} variants")
        print(f"  After:  {after} variants")
        print(f"  Removed: {removed} variants ({removed/before*100:.1f}%)")
        print(f"  ✓ Filtered VCF saved: {os.path.basename(output_vcf)}\n")

    return stats


def print_summary(stats):
    """Print summary statistics."""
    if not stats:
        return
    total_before = sum(s['before'] for s in stats.values())
    total_after = sum(s['after'] for s in stats.values())
    total_removed = sum(s['removed'] for s in stats.values())

    print(f"{'='*60}")
    print(f"SUMMARY: AF Filtering")
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
        f.write("AF Filtering Statistics\n")
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
        description='Filter VCF files by allele fraction (AF) to remove likely germline variants',
        epilog='Example: python 4a_filter_AF.py --input /path/to/Candidates/SNVs-MS1 --af-threshold 0.3'
    )
    parser.add_argument('--input', required=True,
                        help='Input directory containing VCF files')
    parser.add_argument('--af-threshold', type=float, default=0.3,
                        help='Maximum AF to keep (default: 0.3)')
    parser.add_argument('--dry-run', action='store_true',
                        help='Show what would be processed without filtering')

    args = parser.parse_args()

    if not os.path.isdir(args.input):
        print(f"ERROR: Input directory not found: {args.input}")
        return

    if args.af_threshold < 0 or args.af_threshold > 1:
        print(f"ERROR: AF threshold must be between 0 and 1 (got {args.af_threshold})")
        return

    print(f"\n{'='*60}")
    print(f"Step 4a: Allele Fraction (AF) Filtering")
    print(f"{'='*60}\n")

    stats = process_directory(args.input, args.af_threshold, args.dry_run)

    if not args.dry_run:
        print_summary(stats)
        output_dir = os.path.join(args.input, "AF_filtered")
        write_stats_file(stats, output_dir)
        print("✓ AF filtering complete!")
        print(f"  Filtered VCFs saved in: {output_dir}/\n")
    else:
        print("DRY RUN: No files were modified.\n")


if __name__ == "__main__":
    main()

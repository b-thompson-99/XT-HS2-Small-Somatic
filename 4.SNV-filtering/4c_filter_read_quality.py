#!/usr/bin/env python3
"""
Step 4c: Read Quality Filtering

Filters variants based on read-level quality metrics from Mutect2:
- MMQ (Median Mapping Quality): ≥ 40
- ECNT (Event Count within ±50bp): ≤ 2
- MPOS (Median Distance from Read End): ≥ 5

All three criteria must be met to pass.

Input:
    - Directory containing strand-bias-filtered VCF files from step 4b

Output:
    - Filtered VCFs in: {input_dir}/Read_filtered/
    - filtering_stats.txt

Usage:
    python 4c_filter_read_quality.py --input /path/to/SB_filtered
"""

import os
import argparse
from pathlib import Path


def parse_info_field(info_string):
    """Parse VCF INFO field into dictionary."""
    info_dict = {}
    for item in info_string.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            info_dict[key] = value
        else:
            info_dict[item] = True
    return info_dict


def extract_first_value(value_string):
    """Extract first value from comma-separated string."""
    if "," in value_string:
        return value_string.split(",")[0]
    return value_string


def filter_vcf_by_read_quality(input_vcf, output_vcf, min_mmq, max_ecnt, min_mpos):
    """
    Filter VCF by read quality metrics.

    Returns:
        Tuple of (variants_before, variants_after, variants_removed, filter_counts).
    """
    variants_before = 0
    variants_after = 0
    filter_counts = {'low_mmq': 0, 'high_ecnt': 0, 'low_mpos': 0, 'missing_fields': 0}

    with open(input_vcf, 'r') as infile, open(output_vcf, 'w') as outfile:
        for line in infile:
            if line.startswith("##") or line.startswith("#"):
                outfile.write(line)
                continue

            variants_before += 1
            columns = line.strip().split("\t")
            info_dict = parse_info_field(columns[7])

            try:
                ecnt = int(info_dict.get("ECNT", -1))
                mmq = float(extract_first_value(info_dict.get("MMQ", "-1")))
                mpos = int(extract_first_value(info_dict.get("MPOS", "-1")))
            except (ValueError, KeyError, AttributeError):
                filter_counts['missing_fields'] += 1
                continue

            passes_filters = True

            if mmq < min_mmq and mmq != -1:
                filter_counts['low_mmq'] += 1
                passes_filters = False

            if ecnt > max_ecnt and ecnt != -1:
                filter_counts['high_ecnt'] += 1
                passes_filters = False

            if mpos < min_mpos and mpos != -1:
                filter_counts['low_mpos'] += 1
                passes_filters = False

            if mmq == -1 or ecnt == -1 or mpos == -1:
                filter_counts['missing_fields'] += 1
                passes_filters = False

            if passes_filters:
                outfile.write(line)
                variants_after += 1

    variants_removed = variants_before - variants_after
    return variants_before, variants_after, variants_removed, filter_counts


def find_vcf_files(input_dir):
    """Find all VCF files in directory."""
    vcf_files = []
    for filename in os.listdir(input_dir):
        if filename.endswith(".vcf") or filename.endswith(".vcf.gz"):
            vcf_files.append(os.path.join(input_dir, filename))
    return sorted(vcf_files)


def process_directory(input_dir, min_mmq, max_ecnt, min_mpos, dry_run=False):
    """
    Process all VCF files in a directory.

    Returns:
        Dictionary of statistics per sample.
    """
    output_dir = os.path.join(input_dir, "Read_filtered")

    if not dry_run:
        os.makedirs(output_dir, exist_ok=True)

    vcf_files = find_vcf_files(input_dir)

    if len(vcf_files) == 0:
        print(f"ERROR: No VCF files found in {input_dir}")
        return {}

    print(f"Found {len(vcf_files)} VCF files in: {input_dir}")
    print(f"Filter thresholds:")
    print(f"  MMQ ≥ {min_mmq}")
    print(f"  ECNT ≤ {max_ecnt}")
    print(f"  MPOS ≥ {min_mpos}")
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

        before, after, removed, filter_counts = filter_vcf_by_read_quality(
            vcf_file, output_vcf, min_mmq, max_ecnt, min_mpos
        )

        if vcf_file.endswith('.tmp'):
            os.remove(vcf_file)

        stats[sample_name] = {
            'before': before,
            'after': after,
            'removed': removed,
            'filter_counts': filter_counts,
            'percent_kept': (after / before * 100) if before > 0 else 0
        }

        print(f"  Before: {before} variants")
        print(f"  After:  {after} variants")
        print(f"  Removed: {removed} variants ({removed/before*100:.1f}%)")
        print(f"    - Low MMQ: {filter_counts['low_mmq']}")
        print(f"    - High ECNT: {filter_counts['high_ecnt']}")
        print(f"    - Low MPOS: {filter_counts['low_mpos']}")
        print(f"    - Missing fields: {filter_counts['missing_fields']}")
        print(f"  ✓ Filtered VCF saved: {os.path.basename(output_vcf)}\n")

    return stats


def write_stats_file(stats, output_dir):
    """Write filtering statistics to text file."""
    stats_file = os.path.join(output_dir, "filtering_stats.txt")
    total_before = sum(s['before'] for s in stats.values())
    total_after = sum(s['after'] for s in stats.values())
    total_removed = sum(s['removed'] for s in stats.values())
    total_low_mmq = sum(s['filter_counts']['low_mmq'] for s in stats.values())
    total_high_ecnt = sum(s['filter_counts']['high_ecnt'] for s in stats.values())
    total_low_mpos = sum(s['filter_counts']['low_mpos'] for s in stats.values())
    total_missing = sum(s['filter_counts']['missing_fields'] for s in stats.values())

    with open(stats_file, 'w') as f:
        f.write("="*60 + "\n")
        f.write("Read Quality Filtering Statistics\n")
        f.write("="*60 + "\n\n")
        for sample_name, s in sorted(stats.items()):
            f.write(f"{sample_name}:\n")
            f.write(f"  Input variants:  {s['before']}\n")
            f.write(f"  Output variants: {s['after']}\n")
            f.write(f"  Removed:         {s['removed']} ({s['removed']/s['before']*100:.1f}%)\n")
            f.write(f"    - Low MMQ:     {s['filter_counts']['low_mmq']}\n")
            f.write(f"    - High ECNT:   {s['filter_counts']['high_ecnt']}\n")
            f.write(f"    - Low MPOS:    {s['filter_counts']['low_mpos']}\n")
            f.write(f"    - Missing:     {s['filter_counts']['missing_fields']}\n\n")
        f.write("="*60 + "\n")
        f.write("TOTAL:\n")
        f.write(f"  Samples processed: {len(stats)}\n")
        f.write(f"  Input variants:    {total_before}\n")
        f.write(f"  Output variants:   {total_after}\n")
        f.write(f"  Removed:           {total_removed} ({total_removed/total_before*100:.1f}%)\n")
        f.write(f"\nRemoval breakdown:\n")
        f.write(f"  - Low MMQ:   {total_low_mmq}\n")
        f.write(f"  - High ECNT: {total_high_ecnt}\n")
        f.write(f"  - Low MPOS:  {total_low_mpos}\n")
        f.write(f"  - Missing:   {total_missing}\n")
        f.write("="*60 + "\n")

    print(f"✓ Statistics saved to: {stats_file}")


def print_summary(stats):
    """Print summary statistics."""
    if not stats:
        return
    total_before = sum(s['before'] for s in stats.values())
    total_after = sum(s['after'] for s in stats.values())
    total_removed = sum(s['removed'] for s in stats.values())
    total_low_mmq = sum(s['filter_counts']['low_mmq'] for s in stats.values())
    total_high_ecnt = sum(s['filter_counts']['high_ecnt'] for s in stats.values())
    total_low_mpos = sum(s['filter_counts']['low_mpos'] for s in stats.values())
    total_missing = sum(s['filter_counts']['missing_fields'] for s in stats.values())

    print(f"{'='*60}")
    print(f"SUMMARY: Read Quality Filtering")
    print(f"{'='*60}")
    print(f"Samples processed: {len(stats)}")
    print(f"Total variants before: {total_before}")
    print(f"Total variants after:  {total_after}")
    print(f"Total variants removed: {total_removed} ({total_removed/total_before*100:.1f}%)")
    print(f"\nRemoval breakdown:")
    print(f"  - Low MMQ:     {total_low_mmq}")
    print(f"  - High ECNT:   {total_high_ecnt}")
    print(f"  - Low MPOS:    {total_low_mpos}")
    print(f"  - Missing:     {total_missing}")
    print(f"{'='*60}\n")


def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(
        description='Filter VCF files by read quality metrics (MMQ, ECNT, MPOS)',
        epilog='Example: python 4c_filter_read_quality.py --input /path/to/SB_filtered'
    )
    parser.add_argument('--input', required=True,
                        help='Input directory containing strand-bias-filtered VCF files from step 4b')
    parser.add_argument('--min-mmq', type=float, default=40.0,
                        help='Minimum median mapping quality (default: 40)')
    parser.add_argument('--max-ecnt', type=int, default=2,
                        help='Maximum event count (default: 2)')
    parser.add_argument('--min-mpos', type=int, default=5,
                        help='Minimum median distance from read end (default: 5)')
    parser.add_argument('--dry-run', action='store_true',
                        help='Show what would be processed without filtering')

    args = parser.parse_args()

    if not os.path.isdir(args.input):
        print(f"ERROR: Input directory not found: {args.input}")
        return

    print(f"\n{'='*60}")
    print(f"Step 4c: Read Quality Filtering")
    print(f"{'='*60}\n")

    stats = process_directory(args.input, args.min_mmq, args.max_ecnt, args.min_mpos, args.dry_run)

    if not args.dry_run:
        print_summary(stats)
        output_dir = os.path.join(args.input, "Read_filtered")
        write_stats_file(stats, output_dir)
        print("✓ Read quality filtering complete!")
        print(f"  Filtered VCFs saved in: {output_dir}/\n")
    else:
        print("DRY RUN: No files were modified.\n")


if __name__ == "__main__":
    main()

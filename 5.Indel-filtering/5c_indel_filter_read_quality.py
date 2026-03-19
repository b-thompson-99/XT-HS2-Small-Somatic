#!/usr/bin/env python3
"""
Step 5c: Indel Read Quality Filtering

Filters indels based on read-level quality metrics from Mutect2:
- MMQ (Median Mapping Quality): ≥ 40
- ECNT (Event Count within ±50bp): ≤ 2
- MPOS (Median Distance from Read End): ≥ 5

All three criteria must be met to pass.

Input:
    - SB-filtered indel VCFs from step 5b (SB_filtered/)

Output:
    - Filtered VCFs in: {input_dir}/Read_filtered/
    - filtering_stats.txt

Usage:
    python 5c_indel_filter_read_quality.py --input /path/to/Candidates/Indels/AF_filtered/SB_filtered
"""

import os
import argparse
import gzip


def filter_vcf(input_vcf, output_vcf, min_mmq, max_ecnt, min_mpos):
    """
    Filter VCF by read quality metrics.

    Returns:
        Tuple of (variants_before, variants_after, filter_counts).
    """
    variants_before = 0
    variants_after = 0
    filter_counts = {'low_mmq': 0, 'high_ecnt': 0, 'low_mpos': 0, 'missing_fields': 0}

    opener = gzip.open(input_vcf, 'rt') if input_vcf.endswith('.gz') else open(input_vcf, 'r')

    with opener as infile, open(output_vcf, 'w') as outfile:
        for line in infile:
            if line.startswith("#"):
                outfile.write(line)
                continue

            variants_before += 1
            columns = line.strip().split("\t")
            info_field = columns[7]

            info_dict = {}
            for item in info_field.split(";"):
                if "=" in item:
                    key, val = item.split("=", 1)
                    info_dict[key] = val
                else:
                    info_dict[item] = None

            try:
                ecnt = int(info_dict.get("ECNT", -1))
            except ValueError:
                ecnt = -1

            mmq_str = info_dict.get("MMQ", "-1")
            mmq_values = mmq_str.split(",") if mmq_str else ["-1"]
            try:
                mmq = float(mmq_values[0]) if mmq_values[0] != "-1" else -1
            except ValueError:
                mmq = -1

            mpos_str = info_dict.get("MPOS", "-1")
            mpos_values = mpos_str.split(",") if mpos_str else ["-1"]
            try:
                mpos = int(mpos_values[0]) if mpos_values[0] != "-1" else -1
            except ValueError:
                mpos = -1

            passes_filters = True

            if mmq != -1 and mmq < min_mmq:
                filter_counts['low_mmq'] += 1
                passes_filters = False

            if ecnt != -1 and ecnt > max_ecnt:
                filter_counts['high_ecnt'] += 1
                passes_filters = False

            if mpos != -1 and mpos < min_mpos:
                filter_counts['low_mpos'] += 1
                passes_filters = False

            if mmq == -1 or ecnt == -1 or mpos == -1:
                filter_counts['missing_fields'] += 1
                passes_filters = False

            if passes_filters:
                outfile.write(line)
                variants_after += 1

    variants_removed = variants_before - variants_after
    return variants_before, variants_after, filter_counts


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

        before, after, filter_counts = filter_vcf(vcf_file, output_vcf, min_mmq, max_ecnt, min_mpos)
        removed = before - after

        stats[sample_name] = {'before': before, 'after': after, 'removed': removed, 'filter_counts': filter_counts}

        print(f"  Before: {before} variants")
        print(f"  After:  {after} variants")
        print(f"  Removed: {removed} variants ({removed/before*100:.1f}%)")
        print(f"    - Low MMQ: {filter_counts['low_mmq']}")
        print(f"    - High ECNT: {filter_counts['high_ecnt']}")
        print(f"    - Low MPOS: {filter_counts['low_mpos']}")
        print(f"    - Missing fields: {filter_counts['missing_fields']}")
        print()

    return stats


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
    print(f"SUMMARY: Indel Read Quality Filtering")
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
        f.write("Indel Read Quality Filtering Statistics\n")
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


def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(
        description='Filter indel VCF files by read quality metrics (MMQ, ECNT, MPOS)',
        epilog='Example: python 5c_indel_filter_read_quality.py --input /path/to/Candidates/Indels/AF_filtered/SB_filtered'
    )
    parser.add_argument('--input', required=True,
                        help='Input directory containing SB-filtered indel VCF files')
    parser.add_argument('--min-mmq', type=int, default=40,
                        help='Minimum MMQ threshold (default: 40)')
    parser.add_argument('--max-ecnt', type=int, default=2,
                        help='Maximum ECNT threshold (default: 2)')
    parser.add_argument('--min-mpos', type=int, default=5,
                        help='Minimum MPOS threshold (default: 5)')
    parser.add_argument('--dry-run', action='store_true',
                        help='Show what would be processed without filtering')

    args = parser.parse_args()

    if not os.path.isdir(args.input):
        print(f"ERROR: Input directory not found: {args.input}")
        return

    print(f"\n{'='*60}")
    print(f"Step 5c: Indel Read Quality Filtering")
    print(f"{'='*60}\n")

    stats = process_directory(args.input, args.min_mmq, args.max_ecnt, args.min_mpos, args.dry_run)

    if not args.dry_run:
        print_summary(stats)
        output_dir = os.path.join(args.input, "Read_filtered")
        write_stats_file(stats, output_dir)
        print("✓ Indel read quality filtering complete!")
        print(f"  Filtered VCFs saved in: {output_dir}/\n")
    else:
        print("DRY RUN: No files were modified.\n")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Step 4d: Recurrent Artefact Removal

Identifies genomic positions recurrently flagged across control samples and
removes them from patient VCFs. Runs on the login node only — submits a single
SGE job for all compute-intensive work (mpileup across control BAMs).

Method:
1. Generate BED file from all filtered SNV VCFs (MS1 + MS2 + Duplex merged)
2. Run bcftools mpileup on control MS2 BAMs at those positions
3. Identify recurrent loci (positions with alt reads in multiple controls)
4. Filter patient VCFs to remove recurrent positions

Input:
    - Filtered SNV VCFs from step 4c (Read_filtered directories)
    - Control MS2 BAM files and .bai indexes (dedicated directory)

Output:
    - Filtered VCFs in: {input_dir}/Recurrents_removed/
    - Exclusion loci TSV
    - filtering_stats.txt

Requirements:
    - bcftools
    - cyvcf2
    - SGE cluster with qsub

Usage:
    python 4d_recurrent_removal.py \\
        --input-ms1 /path/to/SNVs-MS1/.../Read_filtered \\
        --input-ms2 /path/to/SNVs-MS2/.../Read_filtered \\
        --input-duplex /path/to/SNVs-duplex/.../Read_filtered \\
        --control-bams /path/to/control_MS2_bams \\
        --reference /path/to/reference.fasta

    # Adjust thresholds:
    python 4d_recurrent_removal.py ... --min-alt-reads 3 --min-samples 5

    # Dry run:
    python 4d_recurrent_removal.py ... --dry-run
"""

import os
import sys
import stat
import argparse
import subprocess
from pathlib import Path


# ─── SGE defaults ─────────────────────────────────────────────────────────────
SGE_MEM     = "16G"
SGE_CORES   = 1
SGE_TIME    = "4:00:00"
SGE_JOBNAME = "rec4d"
# ──────────────────────────────────────────────────────────────────────────────


def find_vcf_files(directory):
    """Find all VCF files in a directory."""
    vcf_files = []
    if not os.path.isdir(directory):
        return vcf_files
    for filename in os.listdir(directory):
        if filename.endswith('.vcf') or filename.endswith('.vcf.gz'):
            vcf_files.append(os.path.join(directory, filename))
    return sorted(vcf_files)


def find_bam_files(directory):
    """Find all BAM files in a directory."""
    bam_files = []
    if not os.path.isdir(directory):
        return bam_files
    for filename in os.listdir(directory):
        if filename.endswith('.bam'):
            bam_files.append(os.path.join(directory, filename))
    return sorted(bam_files)


def build_worker_script(args, input_dirs, control_bams):
    """
    Build the self-contained worker Python script that runs inside the SGE job.
    All logic is embedded so the compute node needs no external files beyond
    standard library + cyvcf2.
    """

    worker = f'''#!/usr/bin/env python3
"""
Auto-generated worker for step 4d (Recurrent Artefact Removal).
Do not run directly — submitted by 4d_recurrent_removal.py via qsub.
"""

import os, sys, subprocess, tempfile, shutil, gzip
from collections import defaultdict

try:
    from cyvcf2 import VCF, Writer
    HAS_CYVCF2 = True
except ImportError:
    HAS_CYVCF2 = False
    print("WARNING: cyvcf2 not available, falling back to manual parsing")

# ── Parameters injected at submission time ─────────────────────────────────────
INPUT_DIRS    = {repr(input_dirs)}
CONTROL_BAMS  = {repr(control_bams)}
REFERENCE     = {repr(args.reference)}
BCFTOOLS      = {repr(args.bcftools)}
MIN_ALT_READS = {args.min_alt_reads}
MIN_SAMPLES   = {args.min_samples}
# ──────────────────────────────────────────────────────────────────────────────


def find_vcf_files(directory):
    vcf_files = []
    if not os.path.isdir(directory):
        return vcf_files
    for fn in os.listdir(directory):
        if fn.endswith('.vcf') or fn.endswith('.vcf.gz'):
            vcf_files.append(os.path.join(directory, fn))
    return sorted(vcf_files)


def merge_vcfs_to_bed(vcf_directories, output_bed):
    print("Stage 1: Generating candidate positions BED file")
    positions = set()

    for vcf_dir in vcf_directories:
        vcf_files = find_vcf_files(vcf_dir)
        print(f"  Found {{len(vcf_files)}} VCFs in {{os.path.basename(vcf_dir)}}")
        for vcf_file in vcf_files:
            if HAS_CYVCF2:
                for variant in VCF(vcf_file):
                    positions.add((variant.CHROM, variant.POS))
            else:
                opener = gzip.open(vcf_file, 'rt') if vcf_file.endswith('.gz') else open(vcf_file, 'r')
                with opener as f:
                    for line in f:
                        if line.startswith('#'):
                            continue
                        cols = line.strip().split('\\t')
                        positions.add((cols[0], int(cols[1])))

    with open(output_bed, 'w') as bed:
        for chrom, pos in sorted(positions):
            bed.write(f"{{chrom}}\\t{{pos - 1}}\\t{{pos}}\\n")

    print(f"  Total unique positions: {{len(positions)}}")
    print(f"  BED file: {{output_bed}}\\n")
    return len(positions)


def run_mpileup_on_controls(control_bams, bed_file, reference, bcftools, output_dir):
    print("Stage 2: Running mpileup on control BAMs")
    print(f"  Processing {{len(control_bams)}} control BAMs...")
    os.makedirs(output_dir, exist_ok=True)
    mpileup_vcfs = []

    for bam_file in control_bams:
        sample_name = os.path.basename(bam_file).replace('.bam', '')
        output_vcf  = os.path.join(output_dir, f"{{sample_name}}.mpileup.vcf")

        cmd = (
            f"{{bcftools}} mpileup -f {{reference}} -R {{bed_file}} "
            f"{{bam_file}} -a FORMAT/AD,FORMAT/DP -d 1000000 -Ou | "
            f"{{bcftools}} call -m -Ov -o {{output_vcf}}"
        )
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

        if result.returncode != 0:
            print(f"  WARNING: mpileup failed for {{sample_name}}: {{result.stderr.strip()}}")
            continue

        mpileup_vcfs.append(output_vcf)

    print(f"  Processed {{len(mpileup_vcfs)}}/{{len(control_bams)}} control BAMs successfully\\n")
    return mpileup_vcfs


def parse_mpileup_vcf_manual(vcf_file):
    """Parse bcftools mpileup|call VCF without cyvcf2. Returns (chrom, pos, ref, alt_count) tuples."""
    results = []
    opener = gzip.open(vcf_file, 'rt') if vcf_file.endswith('.gz') else open(vcf_file, 'r')
    with opener as f:
        for line in f:
            if line.startswith('#'):
                continue
            cols = line.strip().split('\\t')
            if len(cols) < 10:
                continue
            chrom, pos, ref, alt = cols[0], int(cols[1]), cols[3], cols[4]
            if alt in ('.', ''):
                continue
            fmt = cols[8].split(':')
            sample = cols[9].split(':')
            if 'AD' not in fmt:
                continue
            ad_idx = fmt.index('AD')
            if ad_idx >= len(sample):
                continue
            ad = sample[ad_idx].split(',')
            if len(ad) < 2:
                continue
            try:
                alt_count = sum(int(x) for x in ad[1:] if x != '.')
            except ValueError:
                continue
            if alt_count > 0:
                results.append((chrom, pos, ref, alt_count))
    return results


def identify_recurrent_loci(mpileup_vcfs, min_alt_support, min_sample_count, output_tsv):
    print("Stage 3: Identifying recurrent loci")
    print(f"  Thresholds: >= {{min_alt_support}} alt reads in >= {{min_sample_count}} samples")

    alt_support_counter = defaultdict(set)

    for vcf_file in mpileup_vcfs:
        sample_name = os.path.basename(vcf_file).replace('.mpileup.vcf', '')

        if HAS_CYVCF2:
            for variant in VCF(vcf_file):
                try:
                    dp4 = variant.INFO.get("DP4")
                    if dp4 and len(dp4) == 4:
                        alt_count = dp4[2] + dp4[3]
                    else:
                        ad = variant.format('AD')
                        if ad is None:
                            continue
                        alt_count = int(sum(ad[0][1:]))
                except Exception:
                    continue
                if alt_count >= min_alt_support:
                    alt_support_counter[(variant.CHROM, variant.POS, variant.REF)].add(sample_name)
        else:
            for chrom, pos, ref, alt_count in parse_mpileup_vcf_manual(vcf_file):
                if alt_count >= min_alt_support:
                    alt_support_counter[(chrom, pos, ref)].add(sample_name)

    exclude_loci = [
        key for key, samples in alt_support_counter.items()
        if len(samples) >= min_sample_count
    ]

    with open(output_tsv, 'w') as tsv:
        tsv.write("CHROM\\tPOS\\tREF\\tALT\\n")
        for chrom, pos, ref in sorted(exclude_loci):
            tsv.write(f"{{chrom}}\\t{{pos}}\\t{{ref}}\\tANY\\n")

    print(f"  Identified {{len(exclude_loci)}} recurrent loci")
    print(f"  Exclusion list: {{output_tsv}}\\n")
    return len(exclude_loci)


def filter_vcfs_by_exclusion_list(input_dir, exclusion_tsv, output_dir):
    exclude_set = set()
    with open(exclusion_tsv, 'r') as tsv:
        next(tsv)
        for line in tsv:
            cols = line.strip().split('\\t')
            if len(cols) >= 3:
                exclude_set.add((cols[0], int(cols[1]), cols[2]))

    vcf_files = find_vcf_files(input_dir)
    os.makedirs(output_dir, exist_ok=True)
    stats = {{}}

    for vcf_file in vcf_files:
        sample_name = os.path.basename(vcf_file).replace('.vcf.gz', '').replace('.vcf', '')
        output_vcf  = os.path.join(output_dir, f"{{sample_name}}.vcf")
        kept_variants    = []
        removed_variants = []

        if HAS_CYVCF2:
            vcf_reader = VCF(vcf_file)
            vcf_writer = Writer(output_vcf, vcf_reader)
            for variant in vcf_reader:
                key = (variant.CHROM, variant.POS, variant.REF, variant.ALT[0] if variant.ALT else '.')
                if (variant.CHROM, variant.POS, variant.REF) in exclude_set:
                    removed_variants.append(key)
                else:
                    vcf_writer.write_record(variant)
                    kept_variants.append(key)
            vcf_writer.close()
        else:
            opener = gzip.open(vcf_file, 'rt') if vcf_file.endswith('.gz') else open(vcf_file, 'r')
            with opener as infile, open(output_vcf, 'w') as outfile:
                for line in infile:
                    if line.startswith('#'):
                        outfile.write(line)
                        continue
                    cols = line.strip().split('\\t')
                    key  = (cols[0], int(cols[1]), cols[3], cols[4])
                    if (cols[0], int(cols[1]), cols[3]) in exclude_set:
                        removed_variants.append(key)
                    else:
                        outfile.write(line)
                        kept_variants.append(key)

        stats[sample_name] = {{
            'kept': len(kept_variants),
            'removed': len(removed_variants),
            'kept_variants': kept_variants,
            'removed_variants': removed_variants
        }}
    return stats


def write_stats_file(stats, output_file, exclusion_tsv=None):
    import datetime
    total_kept    = sum(s['kept']    for s in stats.values())
    total_removed = sum(s['removed'] for s in stats.values())
    total_before  = total_kept + total_removed

    with open(output_file, 'w') as f:
        f.write("=" * 60 + "\\n")
        f.write("Recurrent Artefact Removal Statistics\\n")
        f.write(f"Run: {{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}}\\n")
        f.write("=" * 60 + "\\n\\n")

        if exclusion_tsv and os.path.exists(exclusion_tsv):
            f.write("Recurrent loci excluded:\\n")
            with open(exclusion_tsv, 'r') as tsv:
                next(tsv)
                for line in tsv:
                    cols = line.strip().split('\\t')
                    if len(cols) >= 3:
                        f.write(f"  {{cols[0]}}\\t{{cols[1]}}\\t{{cols[2]}}\\n")
            f.write("\\n")

        for sample_name, s in sorted(stats.items()):
            before = s['kept'] + s['removed']
            pct    = s['removed'] / before * 100 if before > 0 else 0
            f.write(f"{{sample_name}}:\\n")
            f.write(f"  Input variants:  {{before}}\\n")
            f.write(f"  Output variants: {{s['kept']}}\\n")
            f.write(f"  Removed:         {{s['removed']}} ({{pct:.1f}}%)\\n")
            if s['kept_variants']:
                f.write(f"  Kept:\\n")
                for chrom, pos, ref, alt in s['kept_variants']:
                    f.write(f"    {{chrom}}\\t{{pos}}\\t{{ref}}\\t{{alt}}\\n")
            if s['removed_variants']:
                f.write(f"  Removed:\\n")
                for chrom, pos, ref, alt in s['removed_variants']:
                    f.write(f"    {{chrom}}\\t{{pos}}\\t{{ref}}\\t{{alt}}\\n")
            f.write("\\n")

        f.write("=" * 60 + "\\n")
        f.write("TOTAL:\\n")
        f.write(f"  Samples processed: {{len(stats)}}\\n")
        f.write(f"  Input variants:    {{total_before}}\\n")
        f.write(f"  Output variants:   {{total_kept}}\\n")
        pct_total = total_removed / total_before * 100 if total_before > 0 else 0
        f.write(f"  Removed:           {{total_removed}} ({{pct_total:.1f}}%)\\n")
        f.write("=" * 60 + "\\n")


# ── Main ───────────────────────────────────────────────────────────────────────
print(f"\\n{{'='*60}}")
print("Step 4d: Recurrent Artefact Removal")
print(f"{{'='*60}}\\n")

with tempfile.TemporaryDirectory() as temp_dir:

    bed_file      = os.path.join(temp_dir, "candidate_positions.bed")
    num_positions = merge_vcfs_to_bed(list(INPUT_DIRS.values()), bed_file)

    mpileup_dir  = os.path.join(temp_dir, "mpileup_output")
    mpileup_vcfs = run_mpileup_on_controls(
        CONTROL_BAMS, bed_file, REFERENCE, BCFTOOLS, mpileup_dir
    )
    if not mpileup_vcfs:
        print("ERROR: No mpileup VCFs generated -- aborting")
        sys.exit(1)

    exclusion_tsv = os.path.join(temp_dir, "exclusion_loci.tsv")
    num_recurrent = identify_recurrent_loci(
        mpileup_vcfs, MIN_ALT_READS, MIN_SAMPLES, exclusion_tsv
    )

    print("Stage 4: Filtering patient VCFs")
    for label, input_dir in INPUT_DIRS.items():
        output_dir = os.path.join(input_dir, "Recurrents_removed")
        print(f"  Processing {{label}}...")
        stats = filter_vcfs_by_exclusion_list(input_dir, exclusion_tsv, output_dir)
        shutil.copy(exclusion_tsv, os.path.join(output_dir, "exclusion_loci.tsv"))
        write_stats_file(stats, os.path.join(output_dir, "filtering_stats.txt"), exclusion_tsv=exclusion_tsv)
        total_removed = sum(s['removed'] for s in stats.values())
        print(f"    {{len(stats)}} samples, {{total_removed}} variants removed")

print(f"\\n{{'='*60}}")
print("Recurrent Artefact Removal Complete!")
print(f"{{'='*60}}\\n")
print(f"Identified {{num_recurrent}} recurrent artefact positions")
print("\\nFiltered VCFs saved in:")
for label, input_dir in INPUT_DIRS.items():
    print(f"  {{label}}: {{os.path.join(input_dir, 'Recurrents_removed')}}/")
print(f"\\n{{'='*60}}\\n")
print("Step 4 SNV filtering pipeline complete!")
print("Final filtered candidate SNVs are ready for downstream analysis.")
'''
    return worker


def main():
    parser = argparse.ArgumentParser(
        description='Step 4d: Remove recurrent artefacts — submits SGE job',
        epilog='Example: python 4d_recurrent_removal.py --input-ms1 /path/to/Read_filtered '
               '--control-bams /path/to/controls --reference genome.fasta'
    )

    parser.add_argument('--input-ms1',    help='MS1 Read_filtered directory (optional)')
    parser.add_argument('--input-ms2',    help='MS2 Read_filtered directory (optional)')
    parser.add_argument('--input-duplex', help='Duplex Read_filtered directory (optional)')

    parser.add_argument('--control-bams', required=True,
                        help='Directory containing control MS2 BAMs and .bai indexes')
    parser.add_argument('--reference', required=True,
                        help='Reference genome FASTA (must be indexed with samtools faidx)')

    parser.add_argument('--min-alt-reads', type=int, default=5,
                        help='Minimum alt reads in a control to count (default: 5)')
    parser.add_argument('--min-samples', type=int, default=10,
                        help='Minimum controls with alt reads to exclude position (default: 10)')

    parser.add_argument('--bcftools', default='bcftools',
                        help='Path to bcftools executable (default: bcftools)')

    parser.add_argument('--mem',  default=SGE_MEM,
                        help=f'SGE memory request (default: {SGE_MEM})')
    parser.add_argument('--time', default=SGE_TIME,
                        help=f'SGE wallclock time (default: {SGE_TIME})')

    parser.add_argument('--dry-run', action='store_true',
                        help='Show what would be submitted without executing')

    args = parser.parse_args()

    # ── Validate inputs ────────────────────────────────────────────────────────
    input_dirs = {}
    for label, path in [('MS1', args.input_ms1),
                         ('MS2', args.input_ms2),
                         ('Duplex', args.input_duplex)]:
        if path:
            if os.path.isdir(path):
                input_dirs[label] = path
            else:
                print(f"WARNING: {label} directory not found: {path}")

    if not input_dirs:
        print("ERROR: No valid input directories provided.")
        print("Specify at least one of: --input-ms1, --input-ms2, --input-duplex")
        return 1

    if not os.path.isdir(args.control_bams):
        print(f"ERROR: Control BAM directory not found: {args.control_bams}")
        return 1

    control_bams = find_bam_files(args.control_bams)
    if not control_bams:
        print(f"ERROR: No BAM files found in {args.control_bams}")
        return 1

    if not os.path.exists(args.reference):
        print(f"ERROR: Reference genome not found: {args.reference}")
        return 1
    if not os.path.exists(f"{args.reference}.fai"):
        print(f"ERROR: Reference index not found: {args.reference}.fai")
        print(f"  Create with: samtools faidx {args.reference}")
        return 1

    result = subprocess.run([args.bcftools, '--version'], capture_output=True)
    if result.returncode != 0:
        print(f"ERROR: bcftools not found at: {args.bcftools}")
        return 1

    # ── Summary ────────────────────────────────────────────────────────────────
    print(f"\n{'='*60}")
    print("Step 4d: Recurrent Artefact Removal")
    print(f"{'='*60}\n")

    print("Input directories:")
    for label, path in input_dirs.items():
        vcf_count = len(find_vcf_files(path))
        print(f"  {label}: {path} ({vcf_count} VCFs)")

    print(f"\nControl BAMs : {len(control_bams)} files in {args.control_bams}")
    print(f"Reference    : {args.reference}")
    print(f"Thresholds   : >= {args.min_alt_reads} alt reads in >= {args.min_samples} samples")
    print(f"SGE resources: {args.mem} RAM, {args.time} wallclock\n")

    if args.dry_run:
        print("DRY RUN: No job will be submitted.\n")
        return 0

    # ── Build and submit SGE job ───────────────────────────────────────────────
    anchor_dir = list(input_dirs.values())[0]
    log_dir    = os.path.join(os.path.dirname(anchor_dir), "4d_logs")
    os.makedirs(log_dir, exist_ok=True)

    worker_script = os.path.join(log_dir, "4d_worker.py")
    job_script    = os.path.join(log_dir, "4d_job.sh")

    with open(worker_script, 'w') as f:
        f.write(build_worker_script(args, input_dirs, control_bams))

    job_sh = f"""#!/bin/bash
#$ -N {SGE_JOBNAME}
#$ -l h_rt={args.time}
#$ -l mem={args.mem}
#$ -pe smp {SGE_CORES}
#$ -cwd
#$ -o {log_dir}/4d.log
#$ -e {log_dir}/4d.err

echo "Job started: $(date)"
echo "Running on: $(hostname)"
echo ""

python3 {worker_script}

echo ""
echo "Job finished: $(date)"
"""
    with open(job_script, 'w') as f:
        f.write(job_sh)
    os.chmod(job_script, os.stat(job_script).st_mode | stat.S_IXUSR)

    result = subprocess.run(['qsub', job_script], capture_output=True, text=True)

    if result.returncode != 0:
        print(f"ERROR: qsub failed:\n{result.stderr}")
        return 1

    job_id = result.stdout.strip()
    print(f"Job submitted: {job_id}")
    print(f"\nMonitor with:")
    print(f"  qstat")
    print(f"\nLogs:")
    print(f"  {log_dir}/4d.log")
    print(f"  {log_dir}/4d.err")
    print(f"\nOutput VCFs will appear in:")
    for label, path in input_dirs.items():
        print(f"  {label}: {path}/Recurrents_removed/")
    print()

    return 0


if __name__ == "__main__":
    sys.exit(main())

# Step 1c: BWA-MEM Alignment with MBC Preservation

Aligns QC'd paired-end FASTQ files to a reference genome using BWA-MEM. The `-C` flag is critical — it preserves MBC tags in the output BAM. Requires Steps 1a and 1b output.

---

## Input

```
input_fastq_dir/
├── SAMPLE_001/
│   ├── SAMPLE_001_trimmed_R1_fastp.fastq.gz
│   └── SAMPLE_001_trimmed_R2_fastp.fastq.gz
└── ...
```

Reference genome must be indexed with BWA prior to running:
```bash
bwa index /path/to/genome.fasta
```

---

## Usage

Edit `config_1c_alignment.yaml`:

```yaml
paths:
  input_fastq_dir: "/path/to/fastq/data"
  job_scripts_dir: "/path/to/job_scripts"

reference:
  genome_fasta: "/path/to/genome.fasta"

tools:
  bwa_path: "/path/to/bwa"
  samtools_path: "/path/to/samtools"
```

```bash
# Dry run
python 1c_alignment.py --config config_1c_alignment.yaml --dry-run

# Submit jobs
python 1c_alignment.py --config config_1c_alignment.yaml
```

---

## Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `-C` flag | — | Preserves MBC tags in SAM output — do not remove |
| SGE memory | 30G | Per job |
| SGE cores | 8 | Per job |
| SGE wall time | 4:00:00 | Per job |
| Batch size | 5 | Jobs submitted per batch |

---

## Output

```
SAMPLE_001/
├── SAMPLE_001_sorted.bam       # Coordinate-sorted BAM with MBC tags preserved
├── SAMPLE_001_align.stdout
└── SAMPLE_001_align.stderr
```

---

## Troubleshooting

**BWA index missing** — run `bwa index /path/to/genome.fasta` (~1–2 hours for human genome).

**Low mapping rate (<90%)** — verify reference genome version matches sample (hg19 vs hg38); check FASTQ quality from step 1b.

**Out of memory** — increase `sge.memory` in config.

**MBC tags lost in BAM** — verify `-C` flag is present in the BWA command and that input FASTQs contain MBCs in read names (from step 1a).

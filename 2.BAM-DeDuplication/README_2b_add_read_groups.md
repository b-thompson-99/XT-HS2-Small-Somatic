# Step 2b: Add Read Groups to Consensus BAMs

Adds read group (RG) tags to consensus BAMs from step 2a. RG tags are required by Mutect2. Processes all three stringency levels (MS1, MS2, DUPLEX) automatically.

---

## Input

```
creak_output/
├── MS1/SAMPLE_001_hybrid_MS1.bam
├── MS2/SAMPLE_001_hybrid_MS2.bam
└── DUPLEX/SAMPLE_001_duplex.bam
```

---

## Usage

Edit `config_2b_add_read_groups.yaml`:

```yaml
paths:
  creak_output_dir: "/path/to/creak_output"
  job_scripts_dir: "/path/to/job_scripts"

tools:
  picard: "/path/to/picard"
  samtools: "/path/to/samtools"
  java_home: "/path/to/jdk-17"
```

```bash
# Dry run
python 2b_add_read_groups.py --config config_2b_add_read_groups.yaml --dry-run

# Submit jobs
python 2b_add_read_groups.py --config config_2b_add_read_groups.yaml
```

---

## Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `rgpl` | ILLUMINA | Sequencing platform |
| SGE memory | 16G | Per job |
| SGE cores | 4 | Per job |
| SGE wall time | 1:00:00 | Per job |
| Batch size | 5 | Jobs submitted per batch |

`RGSM` (sample name) is extracted automatically from the BAM filename. Other read group fields (`rgid`, `rglb`, `rgpu`) can be left as defaults unless your downstream analysis requires specific values.

---

## Output

```
creak_output/
├── MS1/With_rg/SAMPLE_001_hybrid_MS1.with_rg.bam
│              SAMPLE_001_hybrid_MS1.with_rg.bam.bai
├── MS2/With_rg/SAMPLE_001_hybrid_MS2.with_rg.bam
└── DUPLEX/With_rg/SAMPLE_001_duplex.with_rg.bam
```

> Use the `With_rg/` directories as input for step 3a — not the parent directories.

---

## Troubleshooting

**Picard not found** — verify `picard` path in config.

**No BAM files found** — check step 2a completed successfully and BAMs exist in stringency directories.

**Java errors** — verify Java 17+ is installed and `java_home` is correct in config.

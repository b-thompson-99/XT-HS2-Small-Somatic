# Step 1b: FASTQ QC (fastp)

Runs fastp QC and adapter trimming on paired-end FASTQ files. Requires Step 1a output.

---

## Input

```
input_fastq_dir/
├── SAMPLE_001/
│   ├── SAMPLE_001_trimmed_R1.fastq.gz
│   └── SAMPLE_001_trimmed_R2.fastq.gz
└── ...
```

---

## Usage

Edit `config_1b_fastq_qc.yaml`:

```yaml
paths:
  input_fastq_dir: "/path/to/fastq/data"
  job_scripts_dir: "/path/to/job_scripts"
```

```bash
# Dry run
python 1b_fastq_qc.py --config config_1b_fastq_qc.yaml --dry-run

# Submit jobs
python 1b_fastq_qc.py --config config_1b_fastq_qc.yaml
```

---

## Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| SGE memory | 30G | Per job |
| SGE cores | 4 | Per job |
| SGE wall time | 2:00:00 | Per job |
| Batch size | 5 | Jobs submitted per batch |

---

## Output

```
SAMPLE_001/
├── SAMPLE_001_trimmed_R1_fastp.fastq.gz
├── SAMPLE_001_trimmed_R2_fastp.fastq.gz
├── SAMPLE_001_fastp_report.html
├── SAMPLE_001_fastp_report.json
├── SAMPLE_001_fastp.stdout
└── SAMPLE_001_fastp.stderr
```

Review HTML reports to assess data quality before proceeding to alignment.

---

## Troubleshooting

**No jobs submitted** — verify FASTQ files match the expected naming pattern and sample directories exist in `input_fastq_dir`.

**Jobs fail immediately** — check SGE stderr logs; verify fastp is in PATH.

**Out of memory** — increase `sge.memory` in config.

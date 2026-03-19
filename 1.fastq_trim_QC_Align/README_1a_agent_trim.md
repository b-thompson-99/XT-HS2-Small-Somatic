# Step 1a: MBC Trimming (AGeNT)

Extracts molecular barcodes (MBCs/UMIs) from raw paired-end FASTQ files using Agilent AGeNT trim. 

---

## Input

```
input_fastq_dir/
├── SAMPLE_001/
│   ├── SAMPLE_001_R1.fastq.gz
│   └── SAMPLE_001_R2.fastq.gz
├── SAMPLE_002/
│   ├── SAMPLE_002_R1.fastq.gz
│   └── SAMPLE_002_R2.fastq.gz
└── ...
```

- Paired-end gzipped FASTQs
- Each sample in its own subdirectory
- Files must end with `_R1.fastq.gz` and `_R2.fastq.gz`

---

## Usage

Edit `config_1a_agent_trim.yaml`:

```yaml
paths:
  input_fastq_dir: "/path/to/raw/fastq"
  job_scripts_dir: "/path/to/job_scripts"

tools:
  agent_path: "/path/to/AGeNT/agent/agent.sh"
  java_home: "/path/to/jdk-17"
```

```bash
# Dry run
python 1a_agent_trim.py --config config_1a_agent_trim.yaml --dry-run

# Submit jobs
python 1a_agent_trim.py --config config_1a_agent_trim.yaml
```

---

## Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `adaptor` | `IlluminaXT` | SureSelect XT HS2 adaptor chemistry |
| `mbc_type` | `xths2` | MBC extraction protocol |
| SGE memory | 30G | Per job |
| SGE cores | 6 | Per job |
| SGE wall time | 2:00:00 | Per job |
| Batch size | 5 | Jobs submitted per batch |

---

## Output

```
SAMPLE_001/
├── SAMPLE_001_trimmed_R1.fastq.gz   # Trimmed R1 with MBC in read name
├── SAMPLE_001_trimmed_R2.fastq.gz   # Trimmed R2 with MBC in read name
├── SAMPLE_001_trim.stdout
└── SAMPLE_001_trim.stderr
```

MBC format in read name:
```
@INSTRUMENT:RUNID:FLOWCELL:LANE:TILE:X:Y:MBC_FORWARD:MBC_REVERSE
```

---

## Troubleshooting

**Java errors** — verify Java 17+ is installed (`java -version`) and `java_home` is correct in config.

**AGeNT not found** — verify `agent_path` points to `agent.sh` and the file is executable (`chmod +x`).

**Out of memory** — increase `sge.memory` in config.

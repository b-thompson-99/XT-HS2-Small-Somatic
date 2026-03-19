# Step 2a: CReaK Molecular Consensus Calling (AGeNT)

Groups reads by molecular barcode (MBC) and generates consensus sequences to suppress PCR and sequencing errors. Requires MBC tags from steps 1a–1c.

---

## Stringency levels

| Level | Mode | Description |
|-------|------|-------------|
| MS1 | HYBRID, min 1 | Maximum sensitivity — all UMI families retained |
| MS2 | HYBRID, min 2 | Balanced — minimum 2 reads per UMI family |
| DUPLEX | DUPLEX | Maximum specificity — forward and reverse strand agreement required |

---

## Input

```
input_bam_dir/
├── SAMPLE_001/
│   └── SAMPLE_001_sorted.bam
├── SAMPLE_002/
│   └── SAMPLE_002_sorted.bam
└── ...
```

BAMs must be coordinate-sorted and contain MBC tags in read names (from steps 1a–1c).

---

## Usage

Edit `config_2a_agent_creak.yaml`:

```yaml
paths:
  input_bam_dir: "/path/to/bam/data"
  output_bam_dir: "/path/to/creak_output"
  job_scripts_dir: "/path/to/job_scripts"

tools:
  agent_path: "/path/to/AGeNT/agent/agent.sh"
  java_home: "/path/to/jdk-17"
```

```bash
# Dry run
python 2a_agent_creak.py --config config_2a_agent_creak.yaml --stringency MS1 MS2 DUPLEX --dry-run

# Submit jobs
python 2a_agent_creak.py --config config_2a_agent_creak.yaml --stringency MS1 MS2 DUPLEX
```

The `--stringency` argument accepts one or more values — run all three at once or specify only what is required.

---

## Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `max_distance` | 2 | Max distance (bp) for grouping reads |
| `min_mapping_quality` | 30 | Minimum MAPQ |
| SGE memory | 10G | Per job |
| SGE cores | 6 | Per job |
| SGE wall time | 4:00:00 | Per job |
| Batch size | 5 | Jobs submitted per batch |

---

## Output

```
creak_output/
├── MS1/SAMPLE_001_hybrid_MS1.bam
├── MS2/SAMPLE_001_hybrid_MS2.bam
└── DUPLEX/SAMPLE_001_duplex.bam
```

Consensus BAMs are substantially smaller than input BAMs — this is expected.

---

## Troubleshooting

**Java errors** — verify Java 17+ is installed (`java -version`) and `java_home` is correct in config.

**Out of memory** — increase `sge.memory` in config (try 20–30G); ensure `memory_efficient: true`.

**Very low read count in DUPLEX** — expected for low-coverage samples; consider using MS2 instead. Verify MBC tags are present in input BAMs.

**Empty output BAM** — verify input BAM went through steps 1a–1c and contains MBC tags in read names.

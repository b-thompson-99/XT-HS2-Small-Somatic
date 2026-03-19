# Step 5: Indel Filtering Pipeline

Filters Duplex indel calls from step 3b across four sequential stages.

| Script | Filter | Removes |
|--------|--------|---------|
| `5a_indel_filter_AF.py` | Allele Fraction (AF ≤ 0.3) | Germline variants |
| `5b_indel_filter_SOR.py` | Strand Odds Ratio (SOR ≤ 3) | Strand-biased artefacts |
| `5c_indel_filter_read_quality.py` | MMQ ≥ 40, ECNT ≤ 2, MPOS ≥ 5 | Low-quality calls |
| `5d_indel_exclude_repeats.py` | RepeatMasker overlap | Repetitive region artefacts |

All steps run locally. No SGE submission.

---

## Prerequisites

- Step 3b complete: `Candidates/Indels/` populated (Duplex only)
- RepeatMasker BED file available (see step 5d below)
- `bedtools` installed

---

## Usage

### All four stages via wrapper

```bash
python 5_filter_indels_wrapper.py \
    --input /path/to/Candidates/Indels \
    --repeat-bed /path/to/hg38_repeats.bed
```

Optional parameters (defaults shown):
```bash
python 5_filter_indels_wrapper.py \
    --input /path/to/Candidates/Indels \
    --repeat-bed /path/to/hg38_repeats.bed \
    --af-threshold 0.3 \
    --sor-threshold 3.0 \
    --min-mmq 40 \
    --max-ecnt 2 \
    --min-mpos 5 \
    --dry-run
```

To run steps individually, call each script directly with `--help` for usage.

### Step 5d: RepeatMasker BED file

The RepeatMasker BED file can be downloaded from the UCSC Table Browser (Assembly: hg38, Group: Repeats, Track: RepeatMasker, Output format: BED).

---

## Output Structure

```
Candidates/Indels/
├── SAMPLE.vcf                                          # step 3b output (unfiltered)
└── AF_filtered/
    └── SAMPLE.vcf                                      # after 5a
        └── SB_filtered/
            └── SAMPLE.vcf                              # after 5b
                └── Read_filtered/
                    └── SAMPLE.vcf                      # after 5c
                        └── Repeats_removed/
                            ├── SAMPLE.vcf              # final
                            └── filtering_stats.txt
```

---

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--af-threshold` | 0.3 | Max AF to keep |
| `--sor-threshold` | 3.0 | Max SOR to keep |
| `--min-mmq` | 40 | Min median mapping quality |
| `--max-ecnt` | 2 | Max event count |
| `--min-mpos` | 5 | Min median distance from read end |

---


```

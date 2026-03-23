# Step 4: SNV Filtering Pipeline

Filters Mutect2 SNV calls across four sequential stages to remove false positives and systematic artefacts.

| Script | Filter | Removes |
|--------|--------|---------|
| `4a_filter_AF.py` | Allele Fraction (AF ≤ 0.3) | Germline variants |
| `4b_filter_strand_bias.py` | Strand Odds Ratio (SOR ≤ 3) | Strand-biased artefacts |
| `4c_filter_read_quality.py` | MMQ ≥ 40, ECNT ≤ 2, MPOS ≥ 5 | Low-quality calls |
| `4d_recurrent_removal.py` | Recurrent positions in controls | Systematic sequencing artefacts |

Steps 4a–4c run locally via the wrapper. Step 4d submits an SGE job.

---

## Prerequisites

- Step 3b complete: `Candidates/SNVs-MS1/`, `SNVs-MS2/`, `SNVs-duplex/` populated
- For 4d: control MS2 BAMs and `.bai` indexes in a dedicated directory, `bcftools` available, `cyvcf2` installed

---

## Usage

### Steps 4a–4c: Run wrapper for each stringency

```bash
python 4_filter_snvs_wrapper.py --input /path/to/Candidates/SNVs-MS1
python 4_filter_snvs_wrapper.py --input /path/to/Candidates/SNVs-MS2
python 4_filter_snvs_wrapper.py --input /path/to/Candidates/SNVs-duplex
```

Optional parameters (defaults shown):
```bash
python 4_filter_snvs_wrapper.py \
    --input /path/to/Candidates/SNVs-MS1 \
    --af-threshold 0.3 \
    --sor-threshold 3.0 \
    --min-mmq 40 \
    --max-ecnt 2 \
    --min-mpos 5 \
    --dry-run
```

To run steps individually or adjust thresholds, call `4a_filter_AF.py`, `4b_filter_strand_bias.py`, or `4c_filter_read_quality.py` directly with `--help` for usage.

### Step 4d: Recurrent artefact removal (SGE)

```bash
python 4d_recurrent_removal.py \
    --input-ms1 /path/to/SNVs-MS1/AF_filtered/SB_filtered/Read_filtered \
    --input-ms2 /path/to/SNVs-MS2/AF_filtered/SB_filtered/Read_filtered \
    --input-duplex /path/to/SNVs-duplex/AF_filtered/SB_filtered/Read_filtered \
    --control-bams /path/to/control_MS2_bams \
    --reference /path/to/GRCh38.fasta \
    --dry-run
```

All `--input-*` arguments are optional — specify only the stringencies required.

Default thresholds: `--min-alt-reads 5`, `--min-samples 10`

SGE resources: 16G RAM, 4hr walltime (single job).

---

## Output Structure

```
Candidates/
└── SNVs-MS1/
    ├── SAMPLE.vcf                                      # step 3b output (unfiltered)
    └── AF_filtered/
        └── SAMPLE.vcf                                  # after 4a
            └── SB_filtered/
                └── SAMPLE.vcf                          # after 4b
                    └── Read_filtered/
                        └── SAMPLE.vcf                  # after 4c
                            └── Recurrents_removed/
                                ├── SAMPLE.vcf          # final
                                ├── exclusion_loci.tsv
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
| `--min-alt-reads` | 5 | Min alt reads in a control to count (4d) |
| `--min-samples` | 10 | Min controls with alt reads to exclude position (4d) |

---

## Notes

- On HPC clusters, tools may not be in PATH on compute nodes. If you encounter "command not found" errors, pass full tool paths explicitly — e.g. for step 4d use --bcftools /Full_path/to/bcftools/executable.
- Input VCFs can be `.vcf` or `.vcf.gz`; outputs are uncompressed `.vcf`
- Run on all three stringency levels (MS1, MS2, Duplex)
- Indels have a separate filtering workflow (see `5.Indel-filtering/`)
- 4d uses MS2 consensus BAMs for controls; ensure BAMs have read groups added (step 2b)
- `filtering_stats.txt` in each output directory lists kept and removed variants per sample with CHROM/POS/REF/ALT, plus the full exclusion loci list (4d only)

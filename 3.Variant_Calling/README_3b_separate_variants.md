# Step 3b: Separate SNVs and Indels

Splits Mutect2 VCFs into separate SNV and indel files for downstream filtering. Submits one SGE job per stringency level. Requires step 3a output.

SNVs are extracted from all three stringency levels.
Indels are extracted from Duplex only, as the dual-strand consensus requirement minimises false positives.
Users wishing to call indels at lower stringency would need to modify `3b_separate_variants.py` accordingly.

---

## Input

```
creak_output/
├── MS1/With_rg/Mutect_output/SAMPLE_001.vcf.gz
├── MS2/With_rg/Mutect_output/SAMPLE_001.vcf.gz
└── DUPLEX/With_rg/Mutect_output/SAMPLE_001.vcf.gz
```

---

## Usage

```bash
# Dry run
python 3b_separate_variants.py --base-dir /path/to/creak_output --dry-run

# Submit jobs
python 3b_separate_variants.py --base-dir /path/to/creak_output
```

---

## Output

```
creak_output/Candidates/
├── SNVs-MS1/SAMPLE_001.vcf
├── SNVs-MS2/SAMPLE_001.vcf
├── SNVs-duplex/SAMPLE_001.vcf
└── Indels/SAMPLE_001.vcf
```

Input VCFs (`.vcf` or `.vcf.gz`) are handled automatically. Outputs are uncompressed `.vcf`.

---

## Troubleshooting

**No Mutect_output directories found** — verify step 3a completed and VCFs exist in `With_rg/Mutect_output/`.

**Unexpected variant counts** — expected pattern is MS1 > MS2 > Duplex. Check job logs in `job_scripts/05b_separate/`.

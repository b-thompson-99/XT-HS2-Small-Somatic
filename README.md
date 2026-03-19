# XT-HS2 Small Somatic Variant Calling Pipeline

![Python](https://img.shields.io/badge/python-3.8+-blue) ![Platform](https://img.shields.io/badge/platform-SGE%20cluster-lightgrey)

End-to-end pipeline for calling low-frequency somatic SNVs and indels from Agilent SureSelect XT HS2 sequencing data. Uses UMI consensus collapse (AGeNT CReaK) and GATK Mutect2 to detect variants at allele fractions as low as 0.1%.

Developed for somatic mosaicism studies in post-mortem human brain tissue. Validated on fresh-frozen brain tissue samples from MSA and control cases.

---

## Input Data

Before running the pipeline, organise your FASTQ files as follows — one subdirectory per sample, named after the sample. **Sample names are derived from directory names and propagated through all downstream files**, so consistent naming is important.

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
- Files must end with `_R1.fastq.gz` and `_R2.fastq.gz`
- SGE/UGE cluster with qsub required
- GRCh38 reference genome (no ALT contigs recommended)
- Panel target regions BED file

---

## Dependencies

| Tool | Use |
|------|-----|
| AGeNT (≥ 3.1.2) | FASTQ trimming + CReaK consensus calling |
| fastp | FASTQ QC |
| BWA-MEM (≥ 0.7.17) | Alignment |
| samtools (≥ 1.10) | BAM handling |
| picard | Read group addition |
| GATK (≥ 4.6.1.0) | Mutect2 variant calling |
| Python (≥ 3.8) | All pipeline scripts |
| bcftools | VCF handling (step 4d) |
| cyvcf2 | VCF parsing (step 4d) |
| bedtools | Indel repeat exclusion (step 5d) |

Most tools can be installed via conda (`conda install -c bioconda <tool>`). AGeNT is Agilent software and must be obtained separately — see the AGeNT documentation.

---

## Pipeline Overview

```
Raw paired-end FASTQs
        │
        ▼
  1a. AGeNT Trim ─── MBC/UMI extraction
        │
        ▼
  1b. fastp QC ────── Quality filtering
        │
        ▼
  1c. BWA-MEM ──────── Alignment to GRCh38 (-C flag preserves MBC tags)
        │
        ▼
  2a. AGeNT CReaK ─── UMI consensus collapse (MS1 / MS2 / Duplex)
        │
        ▼
  2b. Picard ──────── Add read groups (required by Mutect2)
        │
        ▼
  3a. GATK Mutect2 ── Somatic variant calling (tumour-only / single sample, no matched normal)
        │
        ▼
  3b. VCF splitting ─ Separate SNVs and indels
        │
    ┌───┴───────────────────────────┐
    │                               │
    ▼                               ▼
SNVs (MS1 + MS2 + Duplex)     Indels (Duplex only)
    │                               │
  4a. AF filter                  5a. AF filter
  4b. Strand bias filter         5b. Strand bias filter
  4c. Read quality filter        5c. Read quality filter
  4d. Recurrent artefact         5d. RepeatMasker exclusion
      removal (SGE)                  │
    │                               ▼
    ▼                        Final candidate indels
Final candidate SNVs
```

---

## Recommendations

We strongly recommend a thorough manual IGV review of all output candidate variants before final analyses. In our validation, duplex-supported candidates were highly reliable; MS1 and MS2 callsets retained some false positives even post-filtering. Manual inspection is therefore highly recommended if using single-strand consensus data for variant calling.

If duplex coverage is sufficient for your study, we would recommend utilising duplex data exclusively. Our duplex coverage was limited by a low library duplication rate, which resulted in small UMI family sizes and reduced the proportion of read pairs eligible for full duplex consensus collapse. For this type of UMI-based duplex sequencing, a high duplication rate is desirable as it reflects larger UMI families and greater potential for duplex agreement.

---

## Quick Start

```bash
# Steps 1a–1c: FASTQ processing
python 1a_agent_trim.py --config config_1a_agent_trim.yaml
python 1b_fastq_qc.py --config config_1b_fastq_qc.yaml
python 1c_alignment.py --config config_1c_alignment.yaml

# Step 2a: CReaK consensus
python 2a_agent_creak.py --config config_2a_agent_creak.yaml --stringency MS1 MS2 DUPLEX

# Step 2b: Read groups
python 2b_add_read_groups.py --config config_2b_add_read_groups.yaml

# Step 3a: Mutect2 (run for each stringency)
python 3a_mutect2_calling.py --config config_3a_mutect2.yaml --input /path/to/creak_output/MS1/With_rg
python 3a_mutect2_calling.py --config config_3a_mutect2.yaml --input /path/to/creak_output/MS2/With_rg
python 3a_mutect2_calling.py --config config_3a_mutect2.yaml --input /path/to/creak_output/DUPLEX/With_rg

# Step 3b: Separate variants
python 3b_separate_variants.py --base-dir /path/to/creak_output

# Step 4: SNV filtering
python 4_filter_snvs_wrapper.py --input /path/to/Candidates/SNVs-MS1
python 4_filter_snvs_wrapper.py --input /path/to/Candidates/SNVs-MS2
python 4_filter_snvs_wrapper.py --input /path/to/Candidates/SNVs-duplex

# Step 4d: Recurrent artefact removal
python 4d_recurrent_removal.py \
    --input-ms1 /path/to/Candidates/SNVs-MS1/AF_filtered/SB_filtered/Read_filtered \
    --input-ms2 /path/to/Candidates/SNVs-MS2/AF_filtered/SB_filtered/Read_filtered \
    --input-duplex /path/to/Candidates/SNVs-duplex/AF_filtered/SB_filtered/Read_filtered \
    --control-bams /path/to/control_MS2_bams \
    --reference /path/to/genome.fasta

# Step 5: Indel filtering
python 5_filter_indels_wrapper.py \
    --input /path/to/Candidates/Indels \
    --repeat-bed /path/to/hg38_repeats.bed
```

For full usage, see [MASTER_COMMAND_LINE_GUIDE.md](MASTER_COMMAND_LINE_GUIDE.md).

---

## Repository Structure

```
XT-HS2-Small-Somatic/
├── 1.fastq_trim_QC_Align/
├── 2.BAM-DeDuplication/
├── 3.Variant_Calling/
├── 4.SNV-filtering/
├── 5.Indel-filtering/
└── MASTER_COMMAND_LINE_GUIDE.md
```

Each directory contains scripts, configs, and a README. All scripts support `--dry-run`.

---

## Output Structure

```
creak_output/
├── MS1/With_rg/Mutect_output/
├── MS2/With_rg/Mutect_output/
├── DUPLEX/With_rg/Mutect_output/
└── Candidates/
    ├── SNVs-MS1/AF_filtered/SB_filtered/Read_filtered/Recurrents_removed/
    ├── SNVs-MS2/AF_filtered/SB_filtered/Read_filtered/Recurrents_removed/
    ├── SNVs-duplex/AF_filtered/SB_filtered/Read_filtered/Recurrents_removed/
    └── Indels/AF_filtered/SB_filtered/Read_filtered/Repeats_removed/
```

---

## Citation

If you use this pipeline, please cite the associated publication (in submission) and the following tools:

- **AGeNT/CReaK**: Agilent Technologies
- **GATK Mutect2**: Benjamin et al. (2019) bioRxiv. doi:10.1101/861054
- **BWA**: Li H. (2013) arXiv:1303.3997v2
- **samtools**: Danecek et al. (2021) *GigaScience* 10(2)
- **bedtools**: Quinlan AR, Hall IM. (2010) *Bioinformatics* 26(6):841–842

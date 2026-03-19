# XT-HS2 Small Somatic Variant Calling Pipeline — Command Line Guide

Complete reference for running the pipeline from raw FASTQs to final somatic SNV and indel candidate callsets.

---

## Full Pipeline Example

```bash
# Steps 1a–1c: FASTQ processing (adapter trimming; QC metrics; alignment to reference genome)
python 1a_agent_trim.py --config config_1a_agent_trim.yaml
python 1b_fastq_qc.py --config config_1b_fastq_qc.yaml
python 1c_alignment.py --config config_1c_alignment.yaml

# Step 2a: CReaK consensus (collapse UMI families into consensus reads at three stringency levels)
python 2a_agent_creak.py --config config_2a_agent_creak.yaml --stringency MS1 MS2 DUPLEX

# Step 2b: Read groups (add read group tags to BAM header required by Mutect2)
python 2b_add_read_groups.py --config config_2b_add_read_groups.yaml

# Step 3a: Mutect2 (somatic SNV and indel calling in tumour-only mode / single sample, no matched normal)
python 3a_mutect2_calling.py --config config_3a_mutect2.yaml --input /data/creak_output/MS1/With_rg
python 3a_mutect2_calling.py --config config_3a_mutect2.yaml --input /data/creak_output/MS2/With_rg
python 3a_mutect2_calling.py --config config_3a_mutect2.yaml --input /data/creak_output/DUPLEX/With_rg

# Step 3b: Separate variants (split Mutect2 VCFs into SNV and indel files)
python 3b_separate_variants.py --base-dir /data/creak_output

# Step 4: SNV filtering (remove germline, strand-biased and low-quality calls)
python 4_filter_snvs_wrapper.py --input /data/creak_output/Candidates/SNVs-MS1
python 4_filter_snvs_wrapper.py --input /data/creak_output/Candidates/SNVs-MS2
python 4_filter_snvs_wrapper.py --input /data/creak_output/Candidates/SNVs-duplex

# Step 4d: Recurrent artefact removal (exclude positions recurrently called in control samples)
python 4d_recurrent_removal.py \
    --input-ms1 /data/creak_output/Candidates/SNVs-MS1/AF_filtered/SB_filtered/Read_filtered \
    --input-ms2 /data/creak_output/Candidates/SNVs-MS2/AF_filtered/SB_filtered/Read_filtered \
    --input-duplex /data/creak_output/Candidates/SNVs-duplex/AF_filtered/SB_filtered/Read_filtered \
    --control-bams /data/control_MS2_bams \
    --reference /data/genome.fasta

# Step 5: Indel filtering (remove germline, strand-biased, low-quality and repeat-region indels)
python 5_filter_indels_wrapper.py \
    --input /data/creak_output/Candidates/Indels \
    --repeat-bed /data/hg38_repeats.bed
```

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
| bcftools | VCF handling |
| cyvcf2 | VCF parsing |
| bedtools | Indel repeat exclusion |

Most tools can be installed via conda (`conda install -c bioconda <tool>`). AGeNT is Agilent software and must be obtained separately — see the AGeNT documentation.

---

## Notes

- All scripts support `--dry-run` to preview actions without executing.
- Input VCFs can be `.vcf` or `.vcf.gz`; outputs are uncompressed `.vcf`.
- For detailed parameter descriptions, see the README in each step's directory.

---

## Directory Structure

```
XT-HS2-Small-Somatic/
├── 1.fastq_trim_QC_Align/
│   ├── 1a_agent_trim.py
│   ├── 1b_fastq_qc.py
│   ├── 1c_alignment.py
│   └── configs + READMEs
├── 2.BAM-DeDuplication/
│   ├── 2a_agent_creak.py
│   ├── 2b_add_read_groups.py
│   └── configs + READMEs
├── 3.Variant_Calling/
│   ├── 3a_mutect2_calling.py
│   ├── 3b_separate_variants.py
│   └── configs + READMEs
├── 4.SNV-filtering/
│   ├── 4_filter_snvs_wrapper.py
│   ├── 4a_filter_AF.py
│   ├── 4b_filter_strand_bias.py
│   ├── 4c_filter_read_quality.py
│   ├── 4d_recurrent_removal.py
│   └── README_4_SNV_filtering.md
├── 5.Indel-filtering/
│   ├── 5_filter_indels_wrapper.py
│   ├── 5a_indel_filter_AF.py
│   ├── 5b_indel_filter_SOR.py
│   ├── 5c_indel_filter_read_quality.py
│   ├── 5d_indel_exclude_repeats.py
│   └── README_5_indel_filtering.md
└── MASTER_COMMAND_LINE_GUIDE.md
```

---

## Steps 1a–1c: FASTQ Trimming, QC and Alignment

Config-based scripts. Edit the relevant YAML file and run:

```bash
python 1a_agent_trim.py --config config_1a_agent_trim.yaml
python 1b_fastq_qc.py --config config_1b_fastq_qc.yaml
python 1c_alignment.py --config config_1c_alignment.yaml
```

Each config requires the following fields to be set:

| Config | Fields to edit |
|--------|---------------|
| `config_1a_agent_trim.yaml` | `input_dir`, `output_dir`, `agent_path` |
| `config_1b_fastq_qc.yaml` | `input_dir`, `output_dir` |
| `config_1c_alignment.yaml` | `input_dir`, `output_dir`, `reference`, `bwa_path` |

See individual READMEs in `1.fastq_trim_QC_Align/` for full config details.

**Output:** =Sorted, indexed BAMs ready for consensus calling.

---

## Step 2a: CReaK Consensus Calling

Generates consensus BAMs from UMI families at three stringency levels:

| Stringency | Description |
|------------|-------------|
| MS1 | Maximum sensitivity — no single-strand multiplicity minimum |
| MS2 | Balanced — single-strand multiplicity minimum of 2 |
| DUPLEX | Maximum specificity — collapses to full duplex consensus (both strands required) |

Edit config — fields to set:

```yaml
# config_2a_agent_creak.yaml
paths:
  input_bam_dir: "/path/to/aligned/bams"   # output from step 1c
  output_dir: "/path/to/creak_output"       # where consensus BAMs will be written
  reference: "/path/to/genome.fasta"        # must be indexed with samtools faidx
tools:
  agent_path: "/path/to/agent"              # path to AGeNT executable
```

Run:
```bash
python 2a_agent_creak.py --config config_2a_agent_creak.yaml --stringency MS1 MS2 DUPLEX
```

The `--stringency` argument accepts one or more values — run all three at once or specify only what is required, e.g. `--stringency DUPLEX` for Duplex only.

**Output:**
```
creak_output/
├── MS1/SAMPLE_hybrid_MS1.bam
├── MS2/SAMPLE_hybrid_MS2.bam
└── DUPLEX/SAMPLE_duplex.bam
```

---

## Step 2b: Add Read Groups

Adds read group tags required by Mutect2.

Edit config — fields to set:

```yaml
# config_2b_add_read_groups.yaml
paths:
  creak_output_dir: "/path/to/creak_output"   # output_dir from step 2a
tools:
  picard: "/path/to/picard"                    # path to Picard jar or executable
  samtools: "/path/to/samtools"
```

Run:
```bash
python 2b_add_read_groups.py --config config_2b_add_read_groups.yaml
```

**Output:** =BAMs with read groups in `With_rg/` subdirectories:
```
creak_output/
├── MS1/With_rg/SAMPLE_hybrid_MS1.with_rg.bam
├── MS2/With_rg/SAMPLE_hybrid_MS2.with_rg.bam
└── DUPLEX/With_rg/SAMPLE_duplex.with_rg.bam
```

> Use the `With_rg/` directories as input for Step 3a.

---

## Step 3a: Mutect2 Variant Calling

Calls somatic SNVs and indels in tumour-only mode (no matched normal). Run separately for each stringency.

Edit config — fields to set:

```yaml
# config_3a_mutect2.yaml
paths:
  reference: "/path/to/genome.fasta"     # must be indexed with samtools faidx and GATK dict
  target_bed: "/path/to/targets.bed"     # panel target regions BED file
tools:
  gatk_path: "/path/to/gatk"             # path to GATK executable or jar
```

Run:
```bash
python 3a_mutect2_calling.py --config config_3a_mutect2.yaml --input /path/to/creak_output/MS1/With_rg
python 3a_mutect2_calling.py --config config_3a_mutect2.yaml --input /path/to/creak_output/MS2/With_rg
python 3a_mutect2_calling.py --config config_3a_mutect2.yaml --input /path/to/creak_output/DUPLEX/With_rg
```

Run for each stringency required — not all three are necessary if only a subset is being used.

**Output:** =Compressed VCFs in `Mutect_output/` within each stringency directory.

---

## Step 3b: Separate SNVs and Indels

Splits Mutect2 VCFs into separate SNV and indel files.

- **SNVs** — all three stringencies (MS1, MS2, Duplex)
- **Indels** — Duplex only

Indels are extracted from Duplex stringency data only, as the dual-strand consensus requirement minimises false positives at low allele fractions. Users wishing to call indels at lower stringency would need to modify `3b_separate_variants.py` accordingly.

```bash
python 3b_separate_variants.py --base-dir /path/to/creak_output
```

**Output:**
```
creak_output/Candidates/
├── SNVs-MS1/
├── SNVs-MS2/
├── SNVs-duplex/
└── Indels/
```

---

## Step 4: SNV Filtering

Four sequential filters applied to each stringency level.

### Steps 4a–4c: Quality filters (wrapper runs all three sequentially)

```bash
python 4_filter_snvs_wrapper.py --input /path/to/Candidates/SNVs-MS1
python 4_filter_snvs_wrapper.py --input /path/to/Candidates/SNVs-MS2
python 4_filter_snvs_wrapper.py --input /path/to/Candidates/SNVs-duplex
```

| Filter | Script | Criterion |
|--------|--------|-----------|
| AF | `4a_filter_AF.py` | AF ≤ 0.3 |
| Strand bias | `4b_filter_strand_bias.py` | SOR ≤ 3 |
| Read quality | `4c_filter_read_quality.py` | MMQ ≥ 40, ECNT ≤ 2, MPOS ≥ 5 |

To run steps individually or adjust thresholds, see `README_4_SNV_filtering.md`.

### Step 4d: Recurrent artefact removal

Identifies alt allele positions recurrently present across control samples and removes them from patient VCFs. Submits a single SGE job — monitor with `qstat`.

Organise control MS2 BAMs and their corresponding `.bai` index files into a dedicated directory, then run:

```bash
python 4d_recurrent_removal.py \
    --input-ms1 /path/to/Candidates/SNVs-MS1/AF_filtered/SB_filtered/Read_filtered \
    --input-ms2 /path/to/Candidates/SNVs-MS2/AF_filtered/SB_filtered/Read_filtered \
    --input-duplex /path/to/Candidates/SNVs-duplex/AF_filtered/SB_filtered/Read_filtered \
    --control-bams /path/to/control_MS2_bams \
    --reference /path/to/genome.fasta
```

All `--input-*` arguments are optional — specify only the stringencies required.

Default thresholds: `--min-alt-reads 5`, `--min-samples 10`

**Final SNV output:**
```
Candidates/SNVs-MS1/AF_filtered/SB_filtered/Read_filtered/Recurrents_removed/
Candidates/SNVs-MS2/AF_filtered/SB_filtered/Read_filtered/Recurrents_removed/
Candidates/SNVs-duplex/AF_filtered/SB_filtered/Read_filtered/Recurrents_removed/
```

---

## Step 5: Indel Filtering

Four sequential filters applied to Duplex indels only. All steps run locally.

### Steps 5a–5d: All filters (wrapper runs all four sequentially)

```bash
python 5_filter_indels_wrapper.py \
    --input /path/to/Candidates/Indels \
    --repeat-bed /path/to/hg38_repeats.bed
```

| Filter | Script | Criterion |
|--------|--------|-----------|
| AF | `5a_indel_filter_AF.py` | AF ≤ 0.3 |
| Strand bias | `5b_indel_filter_SOR.py` | SOR ≤ 3 |
| Read quality | `5c_indel_filter_read_quality.py` | MMQ ≥ 40, ECNT ≤ 2, MPOS ≥ 5 |
| Repeat exclusion | `5d_indel_exclude_repeats.py` | Excludes RepeatMasker regions |

The RepeatMasker BED file (`hg38_repeats.bed`) can be downloaded from the UCSC Table Browser (Assembly: hg38, Group: Repeats, Track: RepeatMasker, Output format: BED).

To run steps individually or adjust thresholds, see `README_5_indel_filtering.md`.

**Final indel output:**
```
Candidates/Indels/AF_filtered/SB_filtered/Read_filtered/Repeats_removed/
```

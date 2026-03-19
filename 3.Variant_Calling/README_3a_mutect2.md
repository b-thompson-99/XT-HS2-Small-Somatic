# Step 3a: Somatic Variant Calling (GATK Mutect2)

Calls somatic SNVs and indels in tumour-only mode (single sample, no matched normal) on consensus BAMs from step 2b. Run separately for each stringency level. Requires step 2b output.

---

## Input

```
creak_output/
├── MS1/With_rg/SAMPLE_001_hybrid_MS1.with_rg.bam
├── MS2/With_rg/SAMPLE_001_hybrid_MS2.with_rg.bam
└── DUPLEX/With_rg/SAMPLE_001_duplex.with_rg.bam
```

Reference genome must be indexed (`samtools faidx genome.fasta`).

---

## Usage

Edit `config_3a_mutect2.yaml`:

```yaml
reference:
  genome_fasta: "/path/to/genome.fasta"
  target_bed: "/path/to/targets.bed"

tools:
  gatk_jar: "/path/to/gatk-package-4.6.1.0-local.jar"
  java_home: "/path/to/jdk-17"
```

```bash
# Dry run
python 3a_mutect2_calling.py --config config_3a_mutect2.yaml --input /path/to/creak_output/MS1/With_rg --dry-run

# Submit jobs — run for each stringency
python 3a_mutect2_calling.py --config config_3a_mutect2.yaml --input /path/to/creak_output/MS1/With_rg
python 3a_mutect2_calling.py --config config_3a_mutect2.yaml --input /path/to/creak_output/MS2/With_rg
python 3a_mutect2_calling.py --config config_3a_mutect2.yaml --input /path/to/creak_output/DUPLEX/With_rg
```

The `--input` flag sets the input directory and automatically creates a `Mutect_output/` subdirectory for outputs.

---

## Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| SGE memory | 9G | Per job — increase for large panels |
| SGE cores | 4 | Per job |
| SGE wall time | 6:00:00 | Per job |

---

## Output

```
With_rg/Mutect_output/
├── SAMPLE_001.vcf.gz
├── SAMPLE_001.vcf.gz.tbi
├── SAMPLE_001_f1r2.tar.gz
├── SAMPLE_001.log
└── SAMPLE_001.err
```

A completion check job is automatically submitted with SGE dependencies — check `check_mutect2_completion.log` to verify all samples completed.

---

## Troubleshooting

**Out of memory** — increase `sge.memory` in config (try 15–20G for large panels).

**Reference index missing** — run `samtools faidx /path/to/genome.fasta`.

**No variants called** — verify the target BED file matches your panel and that BAMs have coverage in target regions.

**GATK version errors** — verify `gatk_jar` path and Java version (`java -version` must be 17+).

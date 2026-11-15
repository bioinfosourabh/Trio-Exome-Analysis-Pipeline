#!/usr/bin/env bash
# Configuration for scripts/germline-variant-calling.sh (trio filters: scripts/trio_filter.sh options)
# Copy this file, fill in the paths for your system, and pass it with -c.
# shellcheck disable=SC2034  # variables are used by the sourcing script

# ---------------------------------------------------------------- reference
# FASTA indexed with: bwa index, samtools faidx, gatk CreateSequenceDictionary
file_ref="/path/to/reference/hg38.fa"

# ---------------------------------------------------------------- known sites
# BQSR known sites. VCFs need an index: .tbi for .vcf.gz, .idx for plain .vcf
known_site="/path/to/resources/Homo_sapiens_assembly38.dbsnp138.vcf"
known_indels=(
    "/path/to/resources/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
    "/path/to/resources/Homo_sapiens_assembly38.known_indels.vcf.gz"
)

# ---------------------------------------------------------------- data
platform="ILLUMINA"
# Target regions (BED) for exome/panel data; leave empty for whole genomes
intervals=""
interval_padding=100

# ---------------------------------------------------------------- filtering
# vqsr: VariantRecalibrator/ApplyVQSR (needs roughly 30+ exomes or WGS)
# hard: GATK generic hard filters (any cohort size)
filter_mode="vqsr"

# VQSR resources (required when filter_mode="vqsr")
training_set_snps="/path/to/resources/hapmap_3.3.hg38.vcf.gz"
training_set_omni="/path/to/resources/1000G_omni2.5.hg38.vcf.gz"            # optional
training_set_1000g="/path/to/resources/1000G_phase1.snps.high_confidence.hg38.vcf.gz"  # optional
training_set_indels="/path/to/resources/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
vqsr_snp_level=99.0
vqsr_indel_level=99.0

# ---------------------------------------------------------------- tools
fastp_args=(--disable_length_filtering --qualified_quality_phred 30)
bwa_args=(-M)
java_mem_small="4g"     # BaseRecalibrator, ApplyBQSR, GenotypeGVCFs
java_mem_large="20g"    # HaplotypeCaller, CombineGVCFs, VQSR
tmp_dir=""              # scratch for GATK/Spark; default: <outdir>/tmp

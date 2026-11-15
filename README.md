# Trio Exome Analysis Pipeline

Germline variant calling followed by trio-based filtering (proband, mother, father) to shortlist candidate variants by inheritance pattern in rare disease and Mendelian disorder studies.

The pipeline has two parts:

1. **Germline variant calling** for all three samples with [`scripts/germiline-variant-calling.sh`](scripts/germiline-variant-calling.sh): fastp, BWA-MEM, GATK BQSR, HaplotypeCaller GVCFs, joint genotyping and VQSR. This is the same script as in [Germline-Variant-Calling-Pipeline](https://github.com/bioinfosourabh/Germline-Variant-Calling-Pipeline); see that repository for its requirements, configuration and outputs.
2. **Trio filtering** with bcftools and awk (commands below), which writes candidate lists for:
   - de novo variants
   - autosomal recessive (homozygous) variants
   - X-linked recessive variants
   - inherited autosomal dominant variants
   - possible parental (maternal) mosaicism

## Requirements

- Everything needed by the germline script (BWA, samtools, fastp, GATK 4.x, bcftools, hg38 reference and resource VCFs)
- `bcftools` and `awk` for the trio filters

## Part 1: Germline variant calling

Put the FASTQ files of the trio (`<sample>_1.fastq.gz`, `<sample>_2.fastq.gz`) in one directory, set the paths in the script (see the Configuration section of the [germline README](https://github.com/bioinfosourabh/Germline-Variant-Calling-Pipeline#configuration)) and run it from that directory:

```sh
cd /path/to/trio_fastq
bash /path/to/scripts/germiline-variant-calling.sh
```

This produces the joint cohort VCF `Germline_VCF/cohort_vqsr.vcf.gz` and per-sample VCFs `<sample>_filtered.vcf.gz`. Functional annotation is not part of the script; if you want annotated per-sample VCFs, annotate them with a tool of your choice (for example VEP or ANNOVAR) and set the file names in Part 2 accordingly.

## Part 2: Trio filtering

The commands below are run in Bash after Part 1.

### Set up sample IDs and paths

```bash
# Sample names as they appear in the VCF header
proband="proband_sample_id"
mother="mother_sample_id"
father="father_sample_id"

# Directory with the per-sample (annotated) VCFs
Dir_vcf="/path/to/vcf/files"

# Per-sample input VCFs (uncompressed; adjust names to your files)
proband_vcf="${Dir_vcf}/${proband}annotated_filterd.vcf"
mother_vcf="${Dir_vcf}/${mother}annotated_filterd.vcf"
father_vcf="${Dir_vcf}/${father}annotated_filterd.vcf"

# Joint VCF containing all three samples
joint_vcf="${Dir_vcf}/cohort_vqsr.vcf.gz"

# Output directory
trio_dir="${Dir_vcf}/Trio_Analysis"
mkdir -p "$trio_dir"
```

### De novo variants

Variants in the proband whose chromosome, position, REF and ALT are absent from both parental VCFs.

```bash
awk 'BEGIN{OFS=FS="\t"} !/^#/ {print $1, $2, $4, $5}' "$mother_vcf" > "$trio_dir/mother_variants.tsv"
awk 'BEGIN{OFS=FS="\t"} !/^#/ {print $1, $2, $4, $5}' "$father_vcf" > "$trio_dir/father_variants.tsv"

awk 'BEGIN{
  FS=OFS="\t"
  while ((getline < "'$trio_dir/mother_variants.tsv'") > 0) mom[$0]=1
  while ((getline < "'$trio_dir/father_variants.tsv'") > 0) dad[$0]=1
}
!/^#/ {
  key = $1"\t"$2"\t"$4"\t"$5
  if (!(key in mom) && !(key in dad)) print
}' "$proband_vcf" > "$trio_dir/${proband}_denovo.vcf"
```

### Autosomal recessive

Homozygous alternate (1/1) in the proband, heterozygous (0/1) in both parents, not on chrX/chrY.

```bash
bcftools view -i 'GT[0]="1/1" && GT[1]="0/1" && GT[2]="0/1" && CHROM !~ "X|Y"' \
    -s "${proband},${mother},${father}" "$joint_vcf" -Oz -o "$trio_dir/${proband}_AR.vcf.gz"
```

### X-linked recessive

On chrX, homozygous or hemizygous alternate in the (male) proband and heterozygous in the mother.

```bash
bcftools view -i 'CHROM="X" && (GT[0]="1/1" || GT[0]="1") && GT[1]="0/1"' \
    -s "${proband},${mother}" "$joint_vcf" -Oz -o "$trio_dir/${proband}_Xlinked.vcf.gz"
```

### Autosomal dominant (inherited)

Heterozygous in the proband and in at least one parent, not on chrX/chrY.

```bash
bcftools view -i 'CHROM !~ "X|Y" && GT[0]="0/1" && (GT[1]="0/1" || GT[2]="0/1")' \
    -s "${proband},${mother},${father}" "$joint_vcf" -Oz -o "$trio_dir/${proband}_Dominant.vcf.gz"
```

### Possible maternal mosaicism

Variants carried by the proband (0/1 or 1/1) where the mother's alternate allele fraction (from AD) is between 10% and 30%. Output columns 5-6 are the proband's GT and AD, and columns 7-8 the mother's.

```bash
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT\t%AD]\n' \
    -s "${proband},${mother}" "$joint_vcf" | \
awk 'BEGIN{OFS="\t"}
{
  split($8, mom_ad, ",");
  mom_ab = mom_ad[2] / (mom_ad[1] + mom_ad[2] + 1e-6);
  if (($5 == "0/1" || $5 == "1/1") && mom_ab > 0.1 && mom_ab < 0.3)
    print $0, "Mosaic_in_mother", mom_ab
}' > "$trio_dir/${proband}_mosaic_in_mother.tsv"

echo "Trio variant filtering completed. Results stored in $trio_dir"
```

## Output

Files in `Trio_Analysis/`:

| File | Content |
|---|---|
| `<proband>_denovo.vcf` | candidate de novo variants |
| `<proband>_AR.vcf.gz` | candidate autosomal recessive variants |
| `<proband>_Xlinked.vcf.gz` | candidate X-linked recessive variants |
| `<proband>_Dominant.vcf.gz` | candidate inherited dominant variants |
| `<proband>_mosaic_in_mother.tsv` | variants with low maternal allele fraction |

## Notes and limitations

- These are genotype-based filters for shortlisting candidates, not a variant classification. Candidates still need annotation, population-frequency filtering and clinical review.
- The genotype filters match unphased genotypes (`0/1`, `1/1`); phased genotypes (`0|1`) are not matched.
- The de novo filter compares sites only, so it does not account for parental sites with no coverage or low genotype quality.
- The mosaicism check is implemented for the mother only; a paternal version needs the father's AD column.
- `GT[0]`, `GT[1]`, `GT[2]` and the awk columns refer to samples by position. Check the sample order of the joint VCF with `bcftools query -l` and make sure it matches proband, mother, father.

## Author

Sourabh Kumar · [bioinfosourabh@gmail.com](mailto:bioinfosourabh@gmail.com) · [bioinfosourabh.netlify.app](https://bioinfosourabh.netlify.app) · ORCID [0009-0009-5840-9830](https://orcid.org/0009-0009-5840-9830)

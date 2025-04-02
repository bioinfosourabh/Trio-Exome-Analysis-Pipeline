# Trio-Exome-Analysis-Pipeline
## 🧬 Germline + Trio Exome Variant Filtering Pipeline
 Author: Sourabh Kumar
This repository contains a step-by-step pipeline to perform germline variant analysis followed by trio-based filtering to identify clinically relevant variant types in rare disease or Mendelian disorder studies.
 Description:
   - Step 1: Perform germline variant calling and annotation for all samples
   - Step 2: Perform trio-based filtering to identify:
       • De novo variants
       • Autosomal Recessive (AR) variants
       • X-linked recessive variants
       • Autosomal dominant variants
       • Mosaic variants in parents


# Germline Analysis Overview
 The initial part of this pipeline performs quality control, alignment,
 BAM processing, variant calling (GVCFs), joint genotyping,
 VQSR (Variant Quality Score Recalibration), and annotation.
 This generates filtered, annotated VCF files for each sample.
 Output: <sample>annotated_filterd.vcf files per sample and joint VCF for cohort.

## Part 0: Set Up Sample Info and File Paths
 Objective:
 Define the proband, mother, and father sample IDs and point to their annotated, filtered VCFs (from previous germline analysis).
Code:
```r
# Define sample names (replace with real sample IDs)
proband="proband_sample_id"
mother="mother_sample_id"
father="father_sample_id"

# Path to directory with filtered/annotated VCFs
Dir_vcf="/path/to/vcf/files"

# Input VCFs (filtered & annotated)
proband_vcf="${Dir_vcf}/${proband}annotated_filterd.vcf"
mother_vcf="${Dir_vcf}/${mother}annotated_filterd.vcf"
father_vcf="${Dir_vcf}/${father}annotated_filterd.vcf"

# Joint VCF with all 3 samples
joint_vcf="${Dir_vcf}/cohort_vqsr.vcf.gz"

# Output directory
trio_dir="${Dir_vcf}/Trio_Analysis"
mkdir -p "$trio_dir"
```

## Part 1: De Novo Variant Detection
Objective:
Find variants that are present in the proband but absent in both parents based on position and alleles.
Code:
```r
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

## Part 2: Autosomal Recessive Inheritance
Objective:
  Find variants that are:
  - Homozygous (1/1) in the proband
  - Heterozygous (0/1) in both parents
  - Located on autosomes (chr1–22)
Code:
```r
bcftools view -i 'GT[0]="1/1" && GT[1]="0/1" && GT[2]="0/1" && CHROM !~ "X|Y"' \
    -s "${proband},${mother},${father}" "$joint_vcf" -Oz -o "$trio_dir/${proband}_AR.vcf.gz"
```

## Part 3: X-Linked Recessive Inheritance
Objective:  
  Find variants that are:
 - On X chromosome
 - Homozygous or hemizygous in male proband
 - Heterozygous in mother
Code:
```r
bcftools view -i 'CHROM="X" && (GT[0]="1/1" || GT[0]="1") && GT[1]="0/1"' \
    -s "${proband},${mother}" "$joint_vcf" -Oz -o "$trio_dir/${proband}_Xlinked.vcf.gz"
```

## Part 4: Autosomal Dominant Inheritance
Objective:
  Find variants that are:
  - Heterozygous (0/1) in proband
  - Also heterozygous in at least one parent
  - Located on autosomes
Code:
```r
bcftools view -i 'CHROM !~ "X|Y" && GT[0]="0/1" && (GT[1]="0/1" || GT[2]="0/1")' \
    -s "${proband},${mother},${father}" "$joint_vcf" -Oz -o "$trio_dir/${proband}_Dominant.vcf.gz"
```

## Part 5: Mosaicism Detection (Parental Mosaic)
Objective:
  Identify cases where:
 - Variant is present in the proband
 - Variant has low allele balance (10–30%) in one parent, suggesting mosaicism

Code:
```r
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%AD]\n' \
    -s "${proband},${mother}" "$joint_vcf" | \
awk 'BEGIN{OFS="\t"}
{
  split($6, mom_ad, ","); 
  mom_ab = mom_ad[2] / (mom_ad[1] + mom_ad[2] + 1e-6);
  if (($5 == "0/1" || $5 == "1/1") && mom_ab > 0.1 && mom_ab < 0.3)
    print $0, "Mosaic_in_mother", mom_ab
}' > "$trio_dir/${proband}_mosaic_in_mother.tsv"

echo "[✔] Trio variant filtering completed. Results stored in $trio_dir"
```

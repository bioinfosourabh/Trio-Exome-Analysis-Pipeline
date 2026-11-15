# Trio Exome Analysis Pipeline

Germline variant calling followed by trio-based filtering (proband, mother, father) to shortlist candidate variants by inheritance pattern in rare disease and Mendelian disorder studies.

The pipeline has two parts:

1. **Germline variant calling** for all three samples with [`scripts/germline-variant-calling.sh`](scripts/germline-variant-calling.sh): fastp, BWA-MEM, GATK BQSR, HaplotypeCaller GVCFs, joint genotyping and VQSR or hard filters. This is the same script as in [Germline-Variant-Calling-Pipeline](https://github.com/bioinfosourabh/Germline-Variant-Calling-Pipeline); see that repository for details of each step.
2. **Trio filtering** of the joint VCF with [`scripts/trio_filter.sh`](scripts/trio_filter.sh) (bcftools and awk), which writes candidate lists for:
   - de novo variants
   - autosomal recessive (homozygous) variants
   - X-linked recessive variants (male proband)
   - inherited autosomal dominant variants
   - possible parental mosaicism (mother and father)

## Requirements

- Everything needed by the germline script (BWA, samtools, fastp, bcftools, tabix, GATK 4.x with Java 17, hg38 reference and resource VCFs)
- `bcftools` and `awk` for the trio filters

## Part 1: Germline variant calling

Put the FASTQ files of the trio (`<sample>_1.fastq.gz`, `<sample>_2.fastq.gz`) in one directory, copy [`config.sh`](config.sh) and set the reference and resource paths in it, then run:

```sh
bash scripts/germline-variant-calling.sh -i /path/to/trio_fastq -o /path/to/results -c my_config.sh -t 16
```

Options: `-i` FASTQ directory, `-o` output directory, `-c` configuration file, `-t` threads, `-s` subset of samples, `-n` dry run, `-f` re-run all steps. Completed steps are skipped when the command is run again, so an interrupted run continues where it stopped. For a trio (three exomes), VQSR has too few variants to train on: set `filter_mode="hard"` in the configuration.

This produces the joint cohort VCF `Germline_VCF/cohort_filtered.vcf.gz`, which Part 2 uses, and per-sample VCFs `Final_VCF/<sample>_filtered.vcf.gz`. Functional annotation is not part of the script; annotate the candidate VCFs from Part 2 with a tool of your choice (for example VEP or ANNOVAR).

## Part 2: Trio filtering

### Set up sample IDs and paths

Sample names are those in the VCF header (`bcftools query -l cohort_filtered.vcf.gz`):

```bash
bash scripts/trio_filter.sh \
    -v /path/to/results/Germline_VCF/cohort_filtered.vcf.gz \
    -p proband_sample_id -m mother_sample_id -f father_sample_id \
    -o /path/to/results/Trio_Analysis \
    -s M            # proband sex: M or F
```

The script checks that the three samples are present, extracts them in the order proband, mother, father (and verifies that order), keeps PASS records (`-a` keeps all), and splits multi-allelic records so that each alternate allele is tested separately. Genotype tests use bcftools genotype keywords, so phased (`0|1`) and unphased (`0/1`) genotypes are treated alike.

### De novo variants

Heterozygous in the proband, homozygous reference in both parents, with at least `-d` reads (default 10) and genotype quality `-g` (default 20) in all three:

```bash
GT[0]="het" && GT[1]="RR" && GT[2]="RR"
  && FMT/DP[0]>=10 && FMT/DP[1]>=10 && FMT/DP[2]>=10
  && FMT/GQ[0]>=20 && FMT/GQ[1]>=20 && FMT/GQ[2]>=20
```

The depth and quality requirements avoid calling a variant de novo because a parent was not covered at that site.

### Autosomal recessive

Homozygous alternate in the proband, heterozygous in both parents, not on chrX, chrY or chrM:

```bash
GT[0]="AA" && GT[1]="het" && GT[2]="het"
```

### X-linked recessive

For a male proband (`-s M`): on chrX, hemizygous or homozygous alternate in the proband, heterozygous in the mother, not carried by the father:

```bash
CHROM="chrX" && (GT[0]="AA" || GT[0]="A") && GT[1]="het" && (GT[2]="RR" || GT[2]="R")
```

### Autosomal dominant (inherited)

Heterozygous in the proband and in at least one parent, autosomes only:

```bash
GT[0]="het" && (GT[1]="het" || GT[2]="het")
```

### Possible parental mosaicism

Variants carried by the proband where a parent's alternate-allele fraction (from `FORMAT/AD`) lies within `-l` (default 0.1 to 0.3), reported separately for the mother and the father.

## Output

Files in the output directory (`-o`):

| File | Content |
|---|---|
| `<proband>_trio.vcf.gz` | trio subset of the joint VCF (proband, mother, father), multi-allelic sites split |
| `<proband>_denovo.vcf.gz` | candidate de novo variants |
| `<proband>_AR.vcf.gz` | candidate autosomal recessive variants |
| `<proband>_Xlinked.vcf.gz` | candidate X-linked recessive variants (male proband) |
| `<proband>_Dominant.vcf.gz` | candidate inherited dominant variants |
| `<proband>_mosaic_in_mother.tsv`, `<proband>_mosaic_in_father.tsv` | variants with a low parental allele fraction |
| `<proband>_trio_summary.tsv` | number of variants per category |

## Notes and limitations

- These are genotype-based filters for shortlisting candidates, not a variant classification. Candidates still need annotation, population-frequency filtering and clinical review.
- Compound heterozygous candidates need gene annotation and are not reported.
- The mosaicism check relies on the allele depths of a germline caller; low-level mosaicism is better assessed with deeper sequencing or a somatic-style caller.
- Version 2.0 replaced the inline commands of earlier versions with `scripts/trio_filter.sh`. Compared with those commands it checks sample names and order instead of assuming them, handles phased genotypes, uses the joint VCF for de novo calls with parental depth and quality checks, and adds paternal mosaicism. The germline script received the fixes described in the [germline repository](https://github.com/bioinfosourabh/Germline-Variant-Calling-Pipeline) (configuration file, argument parsing, input checks, resumable steps, the previously undefined VQSR variables and fastp output folder).

## License

MIT, see [LICENSE](LICENSE).

## Author

Sourabh Kumar · [bioinfosourabh@gmail.com](mailto:bioinfosourabh@gmail.com) · [bioinfosourabh.netlify.app](https://bioinfosourabh.netlify.app) · ORCID [0009-0009-5840-9830](https://orcid.org/0009-0009-5840-9830)

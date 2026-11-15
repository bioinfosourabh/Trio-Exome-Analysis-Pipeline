#!/usr/bin/env bash
#
# Germline short-variant calling (GATK best-practice structure)
#
#   per sample : fastp -> BWA-MEM -> samtools sort -> MarkDuplicatesSpark
#                -> BaseRecalibrator/ApplyBQSR -> HaplotypeCaller (GVCF)
#   cohort     : CombineGVCFs -> GenotypeGVCFs -> VQSR or hard filters
#                -> per-sample VCFs
#
# Each finished step leaves a marker in <outdir>/.steps/; re-running the same
# command skips completed steps and continues where it stopped. A sample that
# fails is reported and left out of joint genotyping; the others continue.
#
# Usage: germline-variant-calling.sh -i FASTQ_DIR -o OUTDIR -c config.sh [-t THREADS]
#
# Author: Sourabh Kumar <bioinfosourabh@gmail.com>
set -euo pipefail

readonly VERSION="2.0.0"
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

usage() {
    cat <<USAGE
Germline variant calling pipeline v${VERSION}

Usage: $(basename "$0") -i FASTQ_DIR -o OUTDIR [-c CONFIG] [options]

  -i DIR     directory with paired FASTQ files <sample>_1.fastq.gz / <sample>_2.fastq.gz  (required)
  -o DIR     output directory                                                           (required)
  -c FILE    configuration file                                     [${script_dir}/../config.sh]
  -t N       threads                                                                    [8]
  -s LIST    comma-separated sample names to run (default: all samples in FASTQ_DIR)
  -n         dry run: check inputs and show the planned steps, run nothing
  -f         force: ignore step markers and re-run everything
  -h         show this help
USAGE
}

input_dir=""; outdir=""; config="${script_dir}/../config.sh"; threads=8; sample_filter=""
dry_run=false; force=false
while getopts ":i:o:c:t:s:nfh" opt; do
    case "$opt" in
        i) input_dir="$OPTARG" ;;
        o) outdir="$OPTARG" ;;
        c) config="$OPTARG" ;;
        t) threads="$OPTARG" ;;
        s) sample_filter="$OPTARG" ;;
        n) dry_run=true ;;
        f) force=true ;;
        h) usage; exit 0 ;;
        :) echo "ERROR: -$OPTARG needs a value" >&2; usage >&2; exit 2 ;;
        *) echo "ERROR: unknown option -$OPTARG" >&2; usage >&2; exit 2 ;;
    esac
done
[[ -n "$input_dir" && -n "$outdir" ]] || { usage >&2; exit 2; }
[[ "$threads" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: -t must be a positive integer" >&2; exit 2; }
[[ -d "$input_dir" ]] || { echo "ERROR: input directory not found: $input_dir" >&2; exit 1; }
[[ -f "$config" ]] || { echo "ERROR: config file not found: $config" >&2; exit 1; }

# ------------------------------------------------------------------ config
# Defaults; the config file overrides them.
file_ref=""; known_site=""; known_indels=(); platform="ILLUMINA"
intervals=""; interval_padding=100; filter_mode="vqsr"
training_set_snps=""; training_set_omni=""; training_set_1000g=""; training_set_indels=""
vqsr_snp_level=99.0; vqsr_indel_level=99.0
fastp_args=(--disable_length_filtering --qualified_quality_phred 30); bwa_args=(-M)
java_mem_small="4g"; java_mem_large="20g"; tmp_dir=""
# shellcheck source=/dev/null
source "$config"

input_dir=$(cd "$input_dir" && pwd)
mkdir -p "$outdir"
outdir=$(cd "$outdir" && pwd)

# Output layout
dir_fastp="$outdir/filtered_qc_report"
dir_map="$outdir/Mapsam"
dir_gvcf="$outdir/GVCFs"
dir_vcf="$outdir/Germline_VCF"
dir_vqsr="$outdir/VQSR"
dir_annotated="$outdir/Final_VCF"
dir_logs="$outdir/logs"
dir_steps="$outdir/.steps"
tmp_dir="${tmp_dir:-$outdir/tmp}"
mkdir -p "$dir_fastp" "$dir_map" "$dir_gvcf" "$dir_vcf" "$dir_vqsr" "$dir_annotated" "$dir_logs" "$dir_steps" "$tmp_dir"

main_log="$dir_logs/pipeline_$(date +%Y%m%d_%H%M%S).log"
log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*" | tee -a "$main_log" >&2; }
die() { log "ERROR: $*"; exit 1; }

# ------------------------------------------------------------------ checks
check_dependencies() {
    local missing=() tool
    for tool in fastp bwa samtools gatk bcftools tabix java; do
        command -v "$tool" >/dev/null 2>&1 || missing+=("$tool")
    done
    if (( ${#missing[@]} > 0 )); then
        if [[ "$dry_run" == true ]]; then
            log "WARNING: not on PATH: ${missing[*]}"
        else
            die "required tools not on PATH: ${missing[*]}"
        fi
    fi
}

vcf_index() {
    if [[ "$1" == *.gz ]]; then echo "$1.tbi"; else echo "$1.idx"; fi
}

require_vcf() {
    local label="$1" path="$2"
    [[ -n "$path" ]] || die "$label is not set in $config"
    [[ -f "$path" ]] || die "$label not found: $path"
    [[ -f "$(vcf_index "$path")" ]] || die "$label has no index ($(vcf_index "$path"))"
}

check_resources() {
    [[ -n "$file_ref" && -f "$file_ref" ]] || die "file_ref not found: '$file_ref'"
    [[ -f "${file_ref}.fai" ]] || die "missing ${file_ref}.fai (samtools faidx)"
    [[ -f "${file_ref%.*}.dict" ]] || die "missing ${file_ref%.*}.dict (gatk CreateSequenceDictionary)"
    [[ -f "${file_ref}.bwt" ]] || die "missing BWA index for $file_ref (bwa index)"
    require_vcf known_site "$known_site"
    local v
    for v in "${known_indels[@]}"; do require_vcf known_indels "$v"; done
    [[ -z "$intervals" || -f "$intervals" ]] || die "intervals file not found: $intervals"
    case "$filter_mode" in
        vqsr)
            require_vcf training_set_snps "$training_set_snps"
            require_vcf training_set_indels "$training_set_indels"
            [[ -z "$training_set_omni" ]] || require_vcf training_set_omni "$training_set_omni"
            [[ -z "$training_set_1000g" ]] || require_vcf training_set_1000g "$training_set_1000g"
            ;;
        hard) ;;
        *) die "filter_mode must be 'vqsr' or 'hard', got '$filter_mode'" ;;
    esac
}

# ------------------------------------------------------------------ step runner
# run_step NAME SCOPE FUNCTION [ARGS...]
# Skips the step if its marker exists and no earlier step of the same scope
# has run in this invocation; otherwise runs FUNCTION with its output in
# logs/SCOPE/NAME.log and writes the marker on success.
declare -A scope_rerun=()
run_step() {
    local name="$1" scope="$2"; shift 2
    local marker="$dir_steps/$scope/$name.done" step_log="$dir_logs/$scope/$name.log"
    mkdir -p "$dir_steps/$scope" "$dir_logs/$scope"
    if [[ -f "$marker" && "$force" == false && -z "${scope_rerun[$scope]:-}" ]]; then
        log "$scope: skip $name (done)"
        return 0
    fi
    if [[ "$dry_run" == true ]]; then
        log "$scope: would run $name"
        return 0
    fi
    scope_rerun[$scope]=1
    rm -f "$marker"
    log "$scope: $name"
    if "$@" > "$step_log" 2>&1; then
        touch "$marker"
    else
        log "$scope: $name FAILED; last lines of $step_log:"
        tail -n 8 "$step_log" | sed 's/^/    /' | tee -a "$main_log" >&2
        return 1
    fi
}

gatk_run() {
    local mem="$1"; shift
    gatk --java-options "-Xmx${mem} -XX:ParallelGCThreads=2 -Djava.io.tmpdir=${tmp_dir}" "$@"
}

interval_args=()
if [[ -n "$intervals" ]]; then
    interval_args=(-L "$intervals" --interval-padding "$interval_padding")
fi

# ------------------------------------------------------------------ per-sample steps
# Step functions chain their commands with && so that a failure anywhere in
# a step is returned to run_step (Bash ignores `set -e` inside functions
# that are called as part of a condition).
step_fastp() {
    local s="$1"
    fastp -i "$input_dir/${s}_1.fastq.gz" -I "$input_dir/${s}_2.fastq.gz" \
        -o "$dir_fastp/${s}_1_filtered.fastq.gz" -O "$dir_fastp/${s}_2_filtered.fastq.gz" \
        --html "$dir_fastp/${s}.html" --json "$dir_fastp/${s}.json" \
        --thread "$(( threads < 16 ? threads : 16 ))" "${fastp_args[@]}"
}

step_align() {
    local s="$1"
    local bam="$dir_map/${s}.bam"
    bwa mem -v 2 -t "$threads" "${bwa_args[@]}" \
        -R "@RG\tID:${s}\tPL:${platform}\tLB:LIB_${s}\tPU:${platform}_${s}\tSM:${s}" \
        "$file_ref" "$dir_fastp/${s}_1_filtered.fastq.gz" "$dir_fastp/${s}_2_filtered.fastq.gz" \
        | samtools sort -@ "$threads" -T "$tmp_dir/${s}.sort" -O BAM -o "$bam" - &&
    samtools index "$bam" &&
    samtools flagstat "$bam" > "$dir_map/${s}.Stat.txt" &&
    if (( $(samtools view -c -F 0x904 "$bam") == 0 )); then
        echo "no mapped reads in $bam" >&2
        return 1
    fi
}

step_markdup() {
    local s="$1"
    gatk_run "$java_mem_large" MarkDuplicatesSpark \
        -I "$dir_map/${s}.bam" -O "$dir_map/${s}_markdup.bam" \
        -M "$dir_map/${s}_markdup_metrics.txt" \
        --spark-master "local[${threads}]" --tmp-dir "$tmp_dir"
}

step_bqsr() {
    local s="$1" ks=(--known-sites "$known_site") v
    for v in "${known_indels[@]}"; do ks+=(--known-sites "$v"); done
    gatk_run "$java_mem_small" BaseRecalibrator -R "$file_ref" \
        -I "$dir_map/${s}_markdup.bam" "${ks[@]}" "${interval_args[@]}" \
        -O "$dir_map/${s}_recall.table" &&
    gatk_run "$java_mem_small" ApplyBQSR -R "$file_ref" \
        -I "$dir_map/${s}_markdup.bam" --bqsr-recal-file "$dir_map/${s}_recall.table" \
        -O "$dir_map/${s}_recall.bam"
}

step_haplotypecaller() {
    local s="$1"
    gatk_run "$java_mem_large" HaplotypeCaller -R "$file_ref" \
        -I "$dir_map/${s}_recall.bam" -O "$dir_gvcf/${s}.g.vcf.gz" -ERC GVCF \
        --native-pair-hmm-threads "$threads" "${interval_args[@]}"
}

process_sample() {
    local s="$1"
    run_step fastp "$s" step_fastp "$s" &&
    run_step align "$s" step_align "$s" &&
    run_step markdup "$s" step_markdup "$s" &&
    run_step bqsr "$s" step_bqsr "$s" &&
    run_step haplotypecaller "$s" step_haplotypecaller "$s"
}

# ------------------------------------------------------------------ cohort steps
step_combine() {
    local vs=() s
    for s in "${done_samples[@]}"; do vs+=(-V "$dir_gvcf/${s}.g.vcf.gz"); done
    gatk_run "$java_mem_large" CombineGVCFs -R "$file_ref" "${vs[@]}" "${interval_args[@]}" \
        -O "$dir_vcf/combined.g.vcf.gz"
}

step_genotype() {
    gatk_run "$java_mem_small" GenotypeGVCFs -R "$file_ref" -V "$dir_vcf/combined.g.vcf.gz" \
        --dbsnp "$known_site" "${interval_args[@]}" -O "$dir_vcf/cohort_raw.vcf.gz"
}

step_vqsr() {
    local snp_res=('--resource:hapmap,known=false,training=true,truth=true,prior=15.0' "$training_set_snps")
    [[ -z "$training_set_omni" ]] || snp_res+=('--resource:omni,known=false,training=true,truth=true,prior=12.0' "$training_set_omni")
    [[ -z "$training_set_1000g" ]] || snp_res+=('--resource:1000G,known=false,training=true,truth=false,prior=10.0' "$training_set_1000g")
    snp_res+=('--resource:dbsnp,known=true,training=false,truth=false,prior=7.0' "$known_site")

    gatk_run "$java_mem_large" VariantRecalibrator -R "$file_ref" -V "$dir_vcf/cohort_raw.vcf.gz" \
        "${snp_res[@]}" -mode SNP -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR \
        -O "$dir_vqsr/recalibrate_SNP.recal" --tranches-file "$dir_vqsr/recalibrate_SNP.tranches" &&
    gatk_run "$java_mem_large" ApplyVQSR -R "$file_ref" -V "$dir_vcf/cohort_raw.vcf.gz" \
        --truth-sensitivity-filter-level "$vqsr_snp_level" -mode SNP \
        --recal-file "$dir_vqsr/recalibrate_SNP.recal" --tranches-file "$dir_vqsr/recalibrate_SNP.tranches" \
        -O "$dir_vqsr/cohort_snps.vcf.gz" &&
    gatk_run "$java_mem_large" VariantRecalibrator -R "$file_ref" -V "$dir_vqsr/cohort_snps.vcf.gz" \
        '--resource:mills,known=false,training=true,truth=true,prior=12.0' "$training_set_indels" \
        '--resource:dbsnp,known=true,training=false,truth=false,prior=2.0' "$known_site" \
        -mode INDEL --max-gaussians 4 -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
        -O "$dir_vqsr/recalibrate_INDEL.recal" --tranches-file "$dir_vqsr/recalibrate_INDEL.tranches" &&
    gatk_run "$java_mem_large" ApplyVQSR -R "$file_ref" -V "$dir_vqsr/cohort_snps.vcf.gz" \
        --truth-sensitivity-filter-level "$vqsr_indel_level" -mode INDEL \
        --recal-file "$dir_vqsr/recalibrate_INDEL.recal" --tranches-file "$dir_vqsr/recalibrate_INDEL.tranches" \
        -O "$dir_vcf/cohort_filtered.vcf.gz"
}

step_hard_filter() {
    local raw="$dir_vcf/cohort_raw.vcf.gz"
    gatk_run "$java_mem_small" SelectVariants -R "$file_ref" -V "$raw" \
        --select-type-to-include SNP -O "$dir_vcf/cohort_snps.vcf.gz" &&
    gatk_run "$java_mem_small" SelectVariants -R "$file_ref" -V "$raw" \
        --select-type-to-include INDEL --select-type-to-include MIXED -O "$dir_vcf/cohort_indels.vcf.gz" &&
    gatk_run "$java_mem_small" VariantFiltration -R "$file_ref" -V "$dir_vcf/cohort_snps.vcf.gz" \
        --filter-name QD2 --filter-expression "QD < 2.0" \
        --filter-name QUAL30 --filter-expression "QUAL < 30.0" \
        --filter-name SOR3 --filter-expression "SOR > 3.0" \
        --filter-name FS60 --filter-expression "FS > 60.0" \
        --filter-name MQ40 --filter-expression "MQ < 40.0" \
        --filter-name MQRankSum-12.5 --filter-expression "MQRankSum < -12.5" \
        --filter-name ReadPosRankSum-8 --filter-expression "ReadPosRankSum < -8.0" \
        -O "$dir_vcf/cohort_snps.filtered.vcf.gz" &&
    gatk_run "$java_mem_small" VariantFiltration -R "$file_ref" -V "$dir_vcf/cohort_indels.vcf.gz" \
        --filter-name QD2 --filter-expression "QD < 2.0" \
        --filter-name QUAL30 --filter-expression "QUAL < 30.0" \
        --filter-name FS200 --filter-expression "FS > 200.0" \
        --filter-name ReadPosRankSum-20 --filter-expression "ReadPosRankSum < -20.0" \
        -O "$dir_vcf/cohort_indels.filtered.vcf.gz" &&
    gatk_run "$java_mem_small" MergeVcfs \
        -I "$dir_vcf/cohort_snps.filtered.vcf.gz" -I "$dir_vcf/cohort_indels.filtered.vcf.gz" \
        -O "$dir_vcf/cohort_filtered.vcf.gz"
}

step_split() {
    local s
    for s in "${done_samples[@]}"; do
        # Sites where the sample carries at least one alternate allele.
        bcftools view -s "$s" --min-ac 1 -Oz -o "$dir_annotated/${s}_filtered.vcf.gz" \
            "$dir_vcf/cohort_filtered.vcf.gz" &&
        tabix -f -p vcf "$dir_annotated/${s}_filtered.vcf.gz" || return 1
    done
}

# ------------------------------------------------------------------ main
log "germline-variant-calling v${VERSION}"
log "input: $input_dir  output: $outdir  config: $config  threads: $threads  filter: $filter_mode"
check_dependencies
check_resources

samples=()
for f in "$input_dir"/*_1.fastq.gz; do
    [[ -e "$f" ]] || continue
    s=$(basename "$f" _1.fastq.gz)
    if [[ ! -f "$input_dir/${s}_2.fastq.gz" ]]; then
        log "WARNING: skipping $s: ${s}_2.fastq.gz not found"
        continue
    fi
    if [[ -n "$sample_filter" && ",${sample_filter}," != *",${s},"* ]]; then
        continue
    fi
    samples+=("$s")
done
(( ${#samples[@]} > 0 )) || die "no FASTQ pairs (<sample>_1.fastq.gz, <sample>_2.fastq.gz) found in $input_dir"
log "samples (${#samples[@]}): ${samples[*]}"

done_samples=(); failed_samples=()
for s in "${samples[@]}"; do
    status=0
    process_sample "$s" || status=$?
    if (( status == 0 )); then
        done_samples+=("$s")
    else
        failed_samples+=("$s")
        log "$s: failed; continuing with the other samples"
    fi
done
(( ${#done_samples[@]} > 0 )) || die "all samples failed"

# A different set of samples invalidates the cohort steps.
cohort_list="$dir_steps/cohort.samples"
if [[ "$dry_run" == false ]]; then
    if [[ -f "$cohort_list" ]] && ! diff -q <(printf '%s\n' "${done_samples[@]}") "$cohort_list" >/dev/null; then
        log "cohort: sample set changed, joint steps will re-run"
        rm -f "$dir_steps/cohort/"*.done
    fi
    printf '%s\n' "${done_samples[@]}" > "$cohort_list"
fi

filter_step=step_hard_filter
[[ "$filter_mode" == vqsr ]] && filter_step=step_vqsr
if ! { run_step combine_gvcfs cohort step_combine &&
        run_step genotype_gvcfs cohort step_genotype &&
        run_step "${filter_mode}_filter" cohort "$filter_step" &&
        run_step split_samples cohort step_split; }; then
    die "joint genotyping failed; see $dir_logs/cohort/"
fi

if [[ "$dry_run" == true ]]; then
    log "dry run complete: nothing was run"
    exit 0
fi
if (( ${#failed_samples[@]} > 0 )); then
    log "finished with failed sample(s): ${failed_samples[*]} (see $dir_logs/<sample>/)"
    exit 1
fi
log "finished: ${#done_samples[@]} sample(s); cohort VCF: $dir_vcf/cohort_filtered.vcf.gz"

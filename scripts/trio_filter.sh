#!/usr/bin/env bash
#
# Trio filtering of a joint-called VCF (proband, mother, father): shortlists
# candidate variants by inheritance pattern.
#
#   de novo               proband het, both parents hom-ref with adequate depth and GQ
#   autosomal recessive   proband hom-alt, both parents het (autosomes)
#   X-linked recessive    male proband hemizygous/hom-alt on chrX, mother het, father not carrying
#   inherited dominant    proband het, at least one parent het (autosomes)
#   parental mosaicism    proband carries the allele, parent's alt fraction within a low range
#
# Genotype tests use bcftools keywords (het, RR, AA, A, R), so phased (0|1)
# and unphased (0/1) genotypes are treated alike. Multi-allelic records are
# split first.
#
# Usage: trio_filter.sh -v cohort.vcf.gz -p PROBAND -m MOTHER -f FATHER -o OUTDIR [options]
#
# Author: Sourabh Kumar <bioinfosourabh@gmail.com>
set -euo pipefail

usage() {
    cat <<USAGE
Usage: $(basename "$0") -v VCF -p PROBAND -m MOTHER -f FATHER -o OUTDIR [options]

  -v FILE    joint-called VCF containing the three samples (bgzipped or plain)  (required)
  -p NAME    proband sample name as in the VCF header                          (required)
  -m NAME    mother sample name                                                (required)
  -f NAME    father sample name                                                (required)
  -o DIR     output directory                                                  (required)
  -s SEX     proband sex, M or F; X-linked recessive is reported for M only    [M]
  -d N       minimum depth (FORMAT/DP) in each trio member for de novo calls    [10]
  -g N       minimum genotype quality (FORMAT/GQ) for de novo calls             [20]
  -l RANGE   parental alt-allele fraction range for mosaicism, LOW,HIGH        [0.1,0.3]
  -a         include records that did not PASS filtering (default: PASS only)
  -h         show this help
USAGE
}

vcf=""; proband=""; mother=""; father=""; outdir=""; sex="M"; min_dp=10; min_gq=20
mosaic="0.1,0.3"; all_records=false
while getopts ":v:p:m:f:o:s:d:g:l:ah" opt; do
    case "$opt" in
        v) vcf="$OPTARG" ;;
        p) proband="$OPTARG" ;;
        m) mother="$OPTARG" ;;
        f) father="$OPTARG" ;;
        o) outdir="$OPTARG" ;;
        s) sex="$OPTARG" ;;
        d) min_dp="$OPTARG" ;;
        g) min_gq="$OPTARG" ;;
        l) mosaic="$OPTARG" ;;
        a) all_records=true ;;
        h) usage; exit 0 ;;
        :) echo "ERROR: -$OPTARG needs a value" >&2; exit 2 ;;
        *) echo "ERROR: unknown option -$OPTARG" >&2; usage >&2; exit 2 ;;
    esac
done
[[ -n "$vcf" && -n "$proband" && -n "$mother" && -n "$father" && -n "$outdir" ]] || { usage >&2; exit 2; }
[[ "$sex" == M || "$sex" == F ]] || { echo "ERROR: -s must be M or F" >&2; exit 2; }
[[ "$min_dp" =~ ^[0-9]+$ && "$min_gq" =~ ^[0-9]+$ ]] || { echo "ERROR: -d and -g must be integers" >&2; exit 2; }
[[ "$mosaic" =~ ^0?\.[0-9]+,0?\.[0-9]+$ ]] || { echo "ERROR: -l must look like 0.1,0.3" >&2; exit 2; }
mosaic_lo=${mosaic%,*}
mosaic_hi=${mosaic#*,}

log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*" >&2; }
die() { log "ERROR: $*"; exit 1; }

command -v bcftools >/dev/null 2>&1 || die "bcftools not found on PATH"
[[ -f "$vcf" ]] || die "VCF not found: $vcf"
mkdir -p "$outdir"

# All three samples must be present, and distinct.
mapfile -t vcf_samples < <(bcftools query -l "$vcf")
for s in "$proband" "$mother" "$father"; do
    printf '%s\n' "${vcf_samples[@]}" | grep -qxF "$s" || die "sample '$s' not in $vcf (samples: ${vcf_samples[*]})"
done
[[ "$proband" != "$mother" && "$proband" != "$father" && "$mother" != "$father" ]] || die "proband, mother and father must be different samples"

# Trio subset in the order proband, mother, father; multi-allelic sites split.
trio="$outdir/${proband}_trio.vcf.gz"
filter_args=()
[[ "$all_records" == true ]] || filter_args=(-f 'PASS,.')
bcftools view -s "${proband},${mother},${father}" "${filter_args[@]}" -Ou "$vcf" \
    | bcftools norm -m -any -Oz -o "$trio"
bcftools index -f -t "$trio"
mapfile -t order < <(bcftools query -l "$trio")
[[ "${order[*]}" == "$proband $mother $father" ]] || die "unexpected sample order in $trio: ${order[*]}"
log "trio VCF: $trio ($(bcftools view -H "$trio" | wc -l) records)"

autosome='CHROM!~"^(chr)?[XYM]T?$"'
chrx='CHROM~"^(chr)?X$"'

select_variants() {
    local name="$1" expr="$2"
    local out="$outdir/${proband}_${name}.vcf.gz"
    bcftools view -i "$expr" -Oz -o "$out" "$trio"
    bcftools index -f -t "$out"
    printf '%s\t%s\n' "$name" "$(bcftools view -H "$out" | wc -l)" >> "$summary"
}

summary="$outdir/${proband}_trio_summary.tsv"
printf 'category\tvariants\n' > "$summary"

select_variants denovo "GT[0]=\"het\" && GT[1]=\"RR\" && GT[2]=\"RR\" && FMT/DP[0]>=${min_dp} && FMT/DP[1]>=${min_dp} && FMT/DP[2]>=${min_dp} && FMT/GQ[0]>=${min_gq} && FMT/GQ[1]>=${min_gq} && FMT/GQ[2]>=${min_gq}"
select_variants AR "${autosome} && GT[0]=\"AA\" && GT[1]=\"het\" && GT[2]=\"het\""
if [[ "$sex" == M ]]; then
    select_variants Xlinked "${chrx} && (GT[0]=\"AA\" || GT[0]=\"A\") && GT[1]=\"het\" && (GT[2]=\"RR\" || GT[2]=\"R\")"
fi
select_variants Dominant "${autosome} && GT[0]=\"het\" && (GT[1]=\"het\" || GT[2]=\"het\")"

# Parental mosaicism: alt-allele fraction of the parent from FORMAT/AD.
mosaic_table() {
    local parent="$1" label="$2"
    local out="$outdir/${proband}_mosaic_in_${label}.tsv"
    {
        printf 'chrom\tpos\tref\talt\tproband_gt\tproband_ad\t%s_gt\t%s_ad\t%s_alt_fraction\n' "$label" "$label" "$label"
        bcftools query -i 'GT[0]="alt"' -f "%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%AD\t]\n" -s "${proband},${parent}" "$trio" \
            | awk -v lo="$mosaic_lo" -v hi="$mosaic_hi" 'BEGIN{OFS="\t"}
                {
                    n = split($8, ad, ",")
                    if (n < 2) next
                    depth = ad[1] + ad[2]
                    if (depth == 0) next
                    frac = ad[2] / depth
                    if (frac > lo && frac < hi) print $1, $2, $3, $4, $5, $6, $7, $8, sprintf("%.3f", frac)
                }'
    } > "$out"
    printf 'mosaic_in_%s\t%s\n' "$label" "$(( $(wc -l < "$out") - 1 ))" >> "$summary"
}
mosaic_table "$mother" mother
mosaic_table "$father" father

log "results in $outdir:"
column -t "$summary" 2>/dev/null >&2 || cat "$summary" >&2

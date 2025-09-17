#!/usr/bin/env bash
# Count per-population callable (n_Loc) and invariant (n_invariant) from an all-sites VCF/BCF.
# Requires: bcftools (>=1.10), GNU parallel
#
# Usage:
#   ./count_sites_per_pop.sh \
#     -v allsites.vcf.gz \
#     -p AFpopfileAutosomal.txt \
#     -o per_population_site_counts.txt \
#     [--to-bcf] [-j 12] [-t 10] [--biallelic-only] [--strict | --max-missing 0.05]
#
# Definitions
#   n_Loc       : callable sites in the population.
#                 default        : AN > 0  (= =1 called genotype in the pop; very relaxed)
#                 --max-missing x: AN = ceil((1 - x) * 2 * n_samples_in_pop)  (e.g., x=0.05 allows =5% missing)
#                 --strict       : AN == 2 * n_samples_in_pop (no missing in the pop; 0% missing)
#   n_invariant : sites where EXACTLY one allele is present among called genotypes
#                 (only REF, or only a single ALT). Multi-allelic handled correctly.

set -euo pipefail

# defaults
JOBS=10
THREADS=8
KEEP_MULTI=1     # if 0 => biallelic-only (-m2 -M2)
STRICT=0
TO_BCF=0
MAX_MISSING="-1" # -1 => disabled; otherwise fraction in [0,1)

VCF=""
POPFILE=""
OUT="per_population_site_counts.txt"

print_help() {
  sed -n '1,120p' "$0" | sed 's/^# \{0,1\}//'
}

# args
while [[ $# -gt 0 ]]; do
  case "$1" in
    -v|--vcf) VCF="$2"; shift 2;;
    -p|--popfile) POPFILE="$2"; shift 2;;
    -o|--out) OUT="$2"; shift 2;;
    -j|--jobs) JOBS="$2"; shift 2;;
    -t|--threads) THREADS="$2"; shift 2;;
    --biallelic-only) KEEP_MULTI=0; shift;;
    --strict) STRICT=1; shift;;
    --max-missing) MAX_MISSING="$2"; shift 2;;
    --to-bcf) TO_BCF=1; shift;;
    -h|--help) print_help; exit 0;;
    *) echo "Unknown arg: $1" >&2; exit 1;;
  esac
done

# checks
command -v bcftools >/dev/null || { echo "bcftools not found" >&2; exit 1; }
command -v parallel  >/dev/null || { echo "GNU parallel not found" >&2; exit 1; }
[[ -n "${VCF}" && -f "${VCF}" ]] || { echo "VCF/BCF missing (-v)"; exit 1; }
[[ -n "${POPFILE}" && -f "${POPFILE}" ]] || { echo "POPFILE missing (-p)"; exit 1; }

# validate missingness/strict combo
if (( STRICT )) && [[ "$MAX_MISSING" != "-1" ]]; then
  echo "Error: --strict and --max-missing are mutually exclusive." >&2
  exit 1
fi
if [[ "$MAX_MISSING" != "-1" ]]; then
  # ensure numeric and in [0,1)
  if ! awk 'BEGIN{v='"$MAX_MISSING"'; if(v>=0 && v<1) exit 0; exit 1}'; then
    echo "Error: --max-missing must be a fraction in [0,1), got '$MAX_MISSING'." >&2
    exit 1
  fi
fi

# normalize popfile line endings
sed -i 's/\r$//' "$POPFILE"
export LC_ALL=C

WORKVCF="$VCF"

# optional: convert to BCF for speed / smaller I/O
if (( TO_BCF )); then
  base="${VCF%.*}"
  outbcf="${base}.bcf"
  if [[ ! -f "$outbcf" ]]; then
    echo "[info] Converting to BCF -> $outbcf"
    bcftools view -Ob -o "$outbcf" "$VCF"
  else
    echo "[info] Using existing BCF: $outbcf"
  fi
  WORKVCF="$outbcf"
fi

# ensure index exists; build if missing
if [[ ! -f "${WORKVCF}.csi" && ! -f "${WORKVCF}.tbi" ]]; then
  echo "[info] Indexing $WORKVCF"
  bcftools index -f "$WORKVCF"
fi

# populations (skip header)
POPS=$(awk -F'\t' 'NR>1{print $2}' "$POPFILE" | sort -u)
[[ -n "$POPS" ]] || { echo "No populations found in $POPFILE" >&2; exit 1; }

# header
printf "Population\tn_Loc\tn_invariant\n" > "$OUT"

export WORKVCF POPFILE THREADS KEEP_MULTI STRICT MAX_MISSING

pop_func() {
  local pop="$1"

  # sample list and count for this pop
  local samples n_samples
  samples=$(awk -F'\t' -v p="$pop" '$2==p{print $1}' "$POPFILE" | paste -sd, -)
  n_samples=$(awk -F'\t' -v p="$pop" '$2==p{c++} END{print c+0}' "$POPFILE")

  if [[ -z "$samples" || "$n_samples" -eq 0 ]]; then
    printf "%s\t0\t0\n" "$pop"
    return
  fi

  # bcftools view (subset to samples; optionally restrict to biallelic)
  local view_cmd
  if (( KEEP_MULTI )); then
    view_cmd=(bcftools view --threads "$THREADS" -s "$samples" -Ou "$WORKVCF")
  else
    view_cmd=(bcftools view --threads "$THREADS" -m2 -M2 -s "$samples" -Ou "$WORKVCF")
  fi

  # view -> fill-tags(AN,AC) -> query -> awk summarize (strict / max-missing / relaxed)
  "${view_cmd[@]}" \
    | bcftools +fill-tags --threads "$THREADS" -Ou -- -t AN,AC \
    | bcftools query -f "%AN\t%AC\n" \
    | awk -F'\t' -v p="$pop" -v strict="$STRICT" -v ns="$n_samples" -v max_missing="$MAX_MISSING" '
        function ceil(x){ return (x==int(x)) ? x : int(x)+1 }
        {
          an = $1 + 0

          if (strict) {
            # strict: no missing in pop
            if (an != 2*ns) next
          } else if (max_missing != "-1") {
            # allow up to max_missing fraction missing
            min_an = (1 - max_missing) * 2 * ns
            min_an = ceil(min_an + 0)   # ensure integer ceiling
            if (an < min_an) next
          } else {
            # relaxed: need at least one call
            if (an <= 0) next
          }

          call++

          # AC can be ".", int, or comma-separated ints (multi-allelic)
          n = split($2, a, ",")
          ac_sum = 0
          alt_present = 0
          for (i=1; i<=n; i++) {
            if (a[i] ~ /^[0-9]+$/) {
              val = a[i] + 0
              ac_sum += val
              if (val > 0) alt_present++
            }
          }
          ref_present = ((an - ac_sum) > 0) ? 1 : 0
          alleles_present = ref_present + alt_present

          # invariant iff exactly one allele present among called genotypes
          if (alleles_present == 1) mono++
        }
        END { printf "%s\t%d\t%d\n", p, call+0, mono+0 }'
}

export -f pop_func

echo "[info] Running $JOBS populations; threads/job=$THREADS; strict=$STRICT; max_missing=${MAX_MISSING}; biallelic_only=$((!KEEP_MULTI))"
printf "%s\n" $POPS | parallel -j "$JOBS" --eta --linebuffer pop_func {} >> "$OUT"

echo "[done] Wrote: $OUT"

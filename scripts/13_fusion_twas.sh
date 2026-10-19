#!/usr/bin/env bash

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/lib/common.sh" "$@"

fusion="$(cfg_path tools.fusion_script)"
ld_prefix="$(cfg_path references.fusion_ld_prefix)"
manifest="$(cfg_path expression.fusion_manifest)"
out="$(output_dir fusion_twas)"
require_file "${fusion}"; require_file "${manifest}"

while IFS=$'\t' read -r tissue position_file weights_dir; do
  [[ "${tissue}" == "tissue" ]] && continue
  position_file="$(resolve_repo_path "${position_file}")"
  weights_dir="$(resolve_repo_path "${weights_dir}")"
  require_file "${position_file}"
  for trait in ra mg; do
    sumstats="$(cfg_path traits.${trait}.file)"
    mkdir -p "${out}/${trait}/${tissue}"
    for chr in $(seq 1 22); do
      Rscript "${fusion}" --sumstats "${sumstats}" --weights "${position_file}" \
        --weights_dir "${weights_dir}" --ref_ld_chr "${ld_prefix}" --chr "${chr}" \
        --out "${out}/${trait}/${tissue}/chr${chr}.dat"
    done
  done
done < "${manifest}"

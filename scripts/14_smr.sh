#!/usr/bin/env bash

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/lib/common.sh" "$@"

smr="$(resolve_tool "$(cfg tools.smr_binary)")"
bfile="$(cfg_path references.smr_bfile_prefix)"
manifest="$(cfg_path expression.smr_manifest)"
out="$(output_dir smr)"
require_file "${manifest}"; require_file "${bfile}.bim"

for trait in ra mg; do
  gwas="${out}/${trait}_gwas.txt"
  Rscript "${REPO_ROOT}/scripts/lib/prepare_smr_gwas.R" "${CONFIG_FILE}" "${trait}" "${gwas}"
  while IFS=$'\t' read -r tissue beqtl_prefix; do
    [[ "${tissue}" == "tissue" ]] && continue
    beqtl_prefix="$(resolve_repo_path "${beqtl_prefix}")"
    "${smr}" --bfile "${bfile}" --gwas-summary "${gwas}" --beqtl-summary "${beqtl_prefix}" \
      --diff-freq-prop 0.1 --out "${out}/${trait}_${tissue}" --thread-num 10
  done < "${manifest}"
done

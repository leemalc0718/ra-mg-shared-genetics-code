#!/usr/bin/env bash

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/lib/common.sh" "$@"

python2="$(resolve_tool "$(cfg tools.python2)")"
hess="$(cfg_path tools.hess_script)"
bfile_pattern="$(cfg_path references.hess_bfile_pattern)"
partition_pattern="$(cfg_path references.hess_partition_pattern)"
out="$(output_dir hess)"
require_file "${hess}"
ra="${out}/ra.hess.tsv.gz"; mg="${out}/mg.hess.tsv.gz"
Rscript "${REPO_ROOT}/scripts/lib/prepare_hess.R" "${CONFIG_FILE}" ra "${ra}"
Rscript "${REPO_ROOT}/scripts/lib/prepare_hess.R" "${CONFIG_FILE}" mg "${mg}"

for chr in $(seq 1 22); do
  bfile="${bfile_pattern//@/${chr}}"
  partition="${partition_pattern//@/${chr}}"
  require_file "${bfile}.bim"; require_file "${partition}"
  "${python2}" "${hess}" --local-rhog "${ra}" "${mg}" --chrom "${chr}" \
    --bfile "${bfile}" --partition "${partition}" --out "${out}/step1"
done

"${python2}" "${hess}" --prefix "${out}/step1_trait1" --out "${out}/step2_trait1"
"${python2}" "${hess}" --prefix "${out}/step1_trait2" --out "${out}/step2_trait2"
"${python2}" "${hess}" --prefix "${out}/step1" \
  --local-hsqg-est "${out}/step2_trait1.txt" "${out}/step2_trait2.txt" \
  --num-shared "$(cfg parameters.hess.number_shared_samples)" \
  --pheno-cor "$(cfg parameters.hess.phenotypic_correlation)" --out "${out}/ra_mg_local_rhog"

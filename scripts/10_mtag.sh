#!/usr/bin/env bash

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/lib/common.sh" "$@"

python2="$(resolve_tool "$(cfg tools.python2)")"
mtag="$(cfg_path tools.mtag_script)"
out="$(output_dir mtag)"
require_file "${mtag}"
Rscript "${REPO_ROOT}/scripts/lib/prepare_mtag.R" "${CONFIG_FILE}" ra "${out}/ra.input.tsv.gz"
Rscript "${REPO_ROOT}/scripts/lib/prepare_mtag.R" "${CONFIG_FILE}" mg "${out}/mg.input.tsv.gz"

"${python2}" "${mtag}" --sumstats "${out}/ra.input.tsv.gz,${out}/mg.input.tsv.gz" \
  --out "${out}/ra_mg" --stream_stdout --use_beta_se --snp_name SNP --a1_name A1 \
  --a2_name A2 --eaf_name eaf --beta_name beta --se_name se --p_name p --n_name n --fdr

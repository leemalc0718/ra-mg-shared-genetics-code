#!/usr/bin/env bash

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/lib/common.sh" "$@"

binary="$(resolve_tool "$(cfg tools.gwas_pw_binary)")"
blocks="$(cfg_path references.gwas_pw_blocks)"
out="$(output_dir gwas_pw)"
input="${out}/ra_mg_input.tsv.gz"
require_file "${blocks}"

Rscript "${REPO_ROOT}/scripts/lib/prepare_gwas_pw.R" "${CONFIG_FILE}" "${input}"
"${binary}" -i "${input}" -bed "${blocks}" -phenos RA MG -cor 0 -o "${out}/ra_mg"

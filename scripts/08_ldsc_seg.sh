#!/usr/bin/env bash

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/lib/common.sh" "$@"

python2="$(resolve_tool "$(cfg tools.python2)")"
ldsc="$(cfg_path tools.ldsc_dir)/ldsc.py"
baseline="$(cfg_path references.ldsc_baseline_prefix)"
weights="$(cfg_path references.ldsc_weights_prefix)"
cts="$(cfg_path references.ldsc_cts_file)"
ldsc_out="$(output_dir ldsc)"
out="$(output_dir ldsc_seg)"
require_file "${ldsc}"; require_file "${cts}"

for trait in ra mg; do
  require_file "${ldsc_out}/${trait}.sumstats.gz"
  "${python2}" "${ldsc}" --h2-cts "${ldsc_out}/${trait}.sumstats.gz" \
    --ref-ld-chr "${baseline}" --ref-ld-chr-cts "${cts}" --w-ld-chr "${weights}" \
    --out "${out}/${trait}_multi_tissue"
done

#!/usr/bin/env bash

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/lib/common.sh" "$@"

python2="$(resolve_tool "$(cfg tools.python2)")"
gnova="$(cfg_path tools.gnova_script)"
bfile="$(cfg_path references.gnova_bfile_pattern)"
ldsc_out="$(output_dir ldsc)"
out="$(output_dir gnova)"
require_file "${gnova}"
require_file "${ldsc_out}/ra.sumstats.gz"
require_file "${ldsc_out}/mg.sumstats.gz"

"${python2}" "${gnova}" "${ldsc_out}/ra.sumstats.gz" "${ldsc_out}/mg.sumstats.gz" \
  --bfile "${bfile}" --out "${out}/ra_mg"

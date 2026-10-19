#!/usr/bin/env bash

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/lib/common.sh" "$@"

python3_bin="$(resolve_tool "$(cfg tools.python3)")"
mixer="$(cfg_path tools.mixer_script)"
figures="$(cfg_path tools.mixer_figures)"
library="$(cfg_path tools.mixer_library)"
bim="$(cfg_path references.mixer_bim_pattern)"
ld="$(cfg_path references.mixer_ld_pattern)"
extract_template="$(cfg_path references.mixer_extract_pattern)"
replicates="$(cfg parameters.mixer.replicates)"
exclude_range="$(cfg parameters.mixer.exclude_range)"
out="$(output_dir mixer)"
require_file "${mixer}"; require_file "${figures}"; require_file "${library}"

Rscript "${REPO_ROOT}/scripts/lib/prepare_mixer.R" "${CONFIG_FILE}" ra "${out}/ra.sumstats.gz"
Rscript "${REPO_ROOT}/scripts/lib/prepare_mixer.R" "${CONFIG_FILE}" mg "${out}/mg.sumstats.gz"

for rep in $(seq 1 "${replicates}"); do
  extract="${extract_template//\{rep\}/${rep}}"
  require_file "${extract}"
  common=(--bim-file "${bim}" --ld-file "${ld}" --extract "${extract}" --lib "${library}" --exclude-ranges "${exclude_range}")
  "${python3_bin}" "${mixer}" fit1 "${common[@]}" --trait1-file "${out}/ra.sumstats.gz" --out "${out}/ra.fit.rep${rep}"
  "${python3_bin}" "${mixer}" fit1 "${common[@]}" --trait1-file "${out}/mg.sumstats.gz" --out "${out}/mg.fit.rep${rep}"
  "${python3_bin}" "${mixer}" fit2 "${common[@]}" --trait1-file "${out}/ra.sumstats.gz" --trait2-file "${out}/mg.sumstats.gz" \
    --trait1-params "${out}/ra.fit.rep${rep}.json" --trait2-params "${out}/mg.fit.rep${rep}.json" --out "${out}/ra_mg.fit.rep${rep}"
  "${python3_bin}" "${mixer}" test1 "${common[@]}" --trait1-file "${out}/ra.sumstats.gz" \
    --load-params "${out}/ra.fit.rep${rep}.json" --out "${out}/ra.test.rep${rep}"
  "${python3_bin}" "${mixer}" test1 "${common[@]}" --trait1-file "${out}/mg.sumstats.gz" \
    --load-params "${out}/mg.fit.rep${rep}.json" --out "${out}/mg.test.rep${rep}"
  "${python3_bin}" "${mixer}" test2 "${common[@]}" --trait1-file "${out}/ra.sumstats.gz" --trait2-file "${out}/mg.sumstats.gz" \
    --load-params "${out}/ra_mg.fit.rep${rep}.json" --out "${out}/ra_mg.test.rep${rep}"
done

"${python3_bin}" "${figures}" combine --json "${out}/ra.fit.rep@.json" --out "${out}/ra.fit"
"${python3_bin}" "${figures}" combine --json "${out}/mg.fit.rep@.json" --out "${out}/mg.fit"
"${python3_bin}" "${figures}" combine --json "${out}/ra_mg.fit.rep@.json" --out "${out}/ra_mg.fit"
"${python3_bin}" "${figures}" combine --json "${out}/ra_mg.test.rep@.json" --out "${out}/ra_mg.test"
"${python3_bin}" "${figures}" one --json "${out}/ra.fit.json" --out "${out}/ra" --statistic mean std
"${python3_bin}" "${figures}" one --json "${out}/mg.fit.json" --out "${out}/mg" --statistic mean std
"${python3_bin}" "${figures}" two --json-fit "${out}/ra_mg.fit.json" --json-test "${out}/ra_mg.test.json" \
  --trait1 RA --trait2 MG --out "${out}/ra_mg" --statistic mean std

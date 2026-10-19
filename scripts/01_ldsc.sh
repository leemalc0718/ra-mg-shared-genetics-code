#!/usr/bin/env bash

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/lib/common.sh" "$@"

python2="$(resolve_tool "$(cfg tools.python2)")"
ldsc_dir="$(cfg_path tools.ldsc_dir)"
munge="${ldsc_dir}/munge_sumstats.py"
ldsc="${ldsc_dir}/ldsc.py"
hm3="$(cfg_path references.hm3_snplist)"
ref_ld="$(cfg_path references.ldsc_eur_prefix)"
ra="$(cfg_path traits.ra.file)"
mg="$(cfg_path traits.mg.file)"
out="$(output_dir ldsc)"

require_file "${munge}"; require_file "${ldsc}"; require_file "${hm3}"
require_file "${ra}"; require_file "${mg}"

"${python2}" "${munge}" --sumstats "${ra}" --snp SNP --a1 A1 --a2 A2 --p p \
  --signed-sumstats beta,0 --frq eaf --N-cas "$(cfg traits.ra.n_case)" \
  --N-con "$(cfg traits.ra.n_control)" --merge-alleles "${hm3}" --out "${out}/ra"

"${python2}" "${munge}" --sumstats "${mg}" --snp SNP --a1 A1 --a2 A2 --p p \
  --signed-sumstats beta,0 --frq eaf --N-cas "$(cfg traits.mg.n_case)" \
  --N-con "$(cfg traits.mg.n_control)" --merge-alleles "${hm3}" --out "${out}/mg"

for trait in ra mg; do
  "${python2}" "${ldsc}" --h2 "${out}/${trait}.sumstats.gz" \
    --ref-ld-chr "${ref_ld}" --w-ld-chr "${ref_ld}" --out "${out}/${trait}_h2"
done

"${python2}" "${ldsc}" --rg "${out}/ra.sumstats.gz,${out}/mg.sumstats.gz" \
  --ref-ld-chr "${ref_ld}" --w-ld-chr "${ref_ld}" --out "${out}/ra_mg_rg"

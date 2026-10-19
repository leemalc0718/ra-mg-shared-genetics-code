#!/usr/bin/env bash

set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

for file in "${ROOT}"/R/*.R "${ROOT}"/scripts/*.R "${ROOT}"/scripts/lib/*.R; do
  Rscript -e "parse(file='${file}')" >/dev/null
done

for file in "${ROOT}"/scripts/*.sh "${ROOT}"/scripts/lib/*.sh; do
  bash -n "${file}"
done

Rscript -e "yaml::read_yaml('${ROOT}/config/config.example.yml')" >/dev/null
python3 "${ROOT}/scripts/lib/config.py" "${ROOT}/config/config.example.yml" traits.ra.n_case >/dev/null

forbidden='/home/|/mnt/|github_pat_|ghp_[A-Za-z0-9]+|Rtmp|wxid_|Users/[A-Za-z0-9._-]+'
if grep -RInE "${forbidden}" "${ROOT}" \
  --exclude='static_checks.sh' --exclude-dir='.git'; then
  echo "Forbidden personal path, temporary path, or credential-like string found." >&2
  exit 1
fi

echo "Static checks passed."

#!/usr/bin/env bash

set -euo pipefail

SCRIPT_LIB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_LIB_DIR}/../.." && pwd)"
CONFIG_READER="${SCRIPT_LIB_DIR}/config.py"

if [[ $# -lt 1 ]]; then
  echo "A configuration file is required." >&2
  exit 2
fi

CONFIG_FILE="$1"
if [[ ! -f "${CONFIG_FILE}" ]]; then
  echo "Configuration file not found: ${CONFIG_FILE}" >&2
  exit 2
fi
CONFIG_FILE="$(cd "$(dirname "${CONFIG_FILE}")" && pwd)/$(basename "${CONFIG_FILE}")"

cfg() {
  python3 "${CONFIG_READER}" "${CONFIG_FILE}" "$1"
}

cfg_path() {
  python3 "${CONFIG_READER}" "${CONFIG_FILE}" "$1" --path
}

resolve_tool() {
  local candidate="$1"
  if [[ "${candidate}" == */* ]]; then
    if [[ "${candidate}" != /* ]]; then candidate="${REPO_ROOT}/${candidate}"; fi
    [[ -e "${candidate}" ]] || { echo "Tool not found: ${candidate}" >&2; exit 2; }
    printf '%s\n' "${candidate}"
  else
    command -v "${candidate}" >/dev/null 2>&1 || { echo "Tool not on PATH: ${candidate}" >&2; exit 2; }
    command -v "${candidate}"
  fi
}

require_file() {
  [[ -f "$1" ]] || { echo "Required file not found: $1" >&2; exit 2; }
}

resolve_repo_path() {
  if [[ "$1" == /* ]]; then printf '%s\n' "$1"; else printf '%s\n' "${REPO_ROOT}/$1"; fi
}

output_dir() {
  local subdir="$1"
  local root
  root="$(cfg_path project.output_dir)"
  mkdir -p "${root}/${subdir}"
  printf '%s\n' "${root}/${subdir}"
}

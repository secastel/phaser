#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_DIR="${ROOT}/data"
TOOLS_DIR="${ROOT}/tools/bin"
if [[ -d "${TOOLS_DIR}" ]]; then
  export PATH="${TOOLS_DIR}:${PATH}"
fi

mkdir -p "${DATA_DIR}"

download() {
  local url="$1"
  local dest="$2"

  if [[ -s "${dest}" ]]; then
    echo "OK: ${dest} (already exists)"
    return 0
  fi

  echo "Downloading: ${url}"
  curl -L --fail --retry 3 --retry-delay 2 -o "${dest}.part" "${url}"
  mv "${dest}.part" "${dest}"
}

# Official upstream test case (large; BAM is multiple GB). URLs from phaser/README.md.
download "https://www.dropbox.com/s/u68p4po2fut2eid/NA06986.vcf.gz?dl=1" "${DATA_DIR}/NA06986.vcf.gz"
download "https://www.dropbox.com/s/328dvei4cqbs7n6/NA06986.vcf.gz.tbi?dl=1" "${DATA_DIR}/NA06986.vcf.gz.tbi"
download "https://www.dropbox.com/s/rxrr01dv4zyhagj/NA06986.2.M_111215_4.bam?dl=1" "${DATA_DIR}/NA06986.2.M_111215_4.bam"
download "https://www.dropbox.com/s/vunmr97j8v6dqi8/NA06986.2.M_111215_4.bam.bai?dl=1" "${DATA_DIR}/NA06986.2.M_111215_4.bam.bai"

test -s "${DATA_DIR}/NA06986.vcf.gz"
test -s "${DATA_DIR}/NA06986.vcf.gz.tbi"
test -s "${DATA_DIR}/NA06986.2.M_111215_4.bam"
test -s "${DATA_DIR}/NA06986.2.M_111215_4.bam.bai"

echo "Done. Test case is in: ${DATA_DIR}"


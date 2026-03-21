#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_DIR="${ROOT}/data"
OUT_DIR="${ROOT}/out"
TOOLS_DIR="${ROOT}/tools/bin"
if [[ -d "${TOOLS_DIR}" ]]; then
  export PATH="${TOOLS_DIR}:${PATH}"
fi

VCF="${DATA_DIR}/NA06986.vcf.gz"
BAM="${DATA_DIR}/NA06986.2.M_111215_4.bam"

if [[ ! -s "${VCF}" || ! -s "${VCF}.tbi" || ! -s "${BAM}" ]]; then
  echo "Missing NA06986 test data in ${DATA_DIR}."
  echo "Run: ./scripts/download_test_case.sh"
  exit 1
fi

mkdir -p "${OUT_DIR}"

require_tools() {
  local missing=0
  for exe in samtools bedtools bcftools bgzip tabix; do
    if ! command -v "${exe}" >/dev/null 2>&1; then
      echo "Missing dependency: ${exe}"
      missing=1
    fi
  done
  if [[ "${missing}" -ne 0 ]]; then
    echo "Install dependencies or add ROOT/tools/bin to PATH."
    exit 1
  fi
}

require_tools

PREFIX="${OUT_DIR}/phaser_test_case"
TIME_FILE="${PREFIX}.time.txt"

TIME_CMD="/usr/bin/time -p"
if [[ ! -x "/usr/bin/time" ]]; then
  TIME_CMD="time"
fi

echo "Running phASER (single thread; default MAPQ=255)..."
echo "Output prefix: ${PREFIX}"

# Keep MAPQ default (255) and single-thread settings for apples-to-apples benchmarking.
${TIME_CMD} python3 "${ROOT}/phaser/phaser.py" \
  --threads 1 --io_threads 0 \
  --vcf "${VCF}" \
  --bam "${BAM}" \
  --paired_end 1 \
  --mapq 255 --baseq 10 \
  --sample NA06986 \
  --prefilter_hets 1 --reintegrate_vcf 1 \
  --o "${PREFIX}" 2> "${TIME_FILE}"

echo "Done."
if [[ -s "${TIME_FILE}" ]] && grep -q "^real" "${TIME_FILE}"; then
  echo "Runtime:"
  cat "${TIME_FILE}"
else
  echo "Runtime saved to: ${TIME_FILE}"
fi


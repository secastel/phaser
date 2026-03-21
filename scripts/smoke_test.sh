#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_DIR="${ROOT}/data"
OUT_DIR="${ROOT}/out/smoke"
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

TIME_CMD="/usr/bin/time -p"
if [[ ! -x "/usr/bin/time" ]]; then
  TIME_CMD="time"
fi

echo "== phASER smoke test (chr1) =="
PHASER_PREFIX="${OUT_DIR}/phaser_chr1"
${TIME_CMD} python3 "${ROOT}/phaser/phaser.py" \
  --threads 1 --io_threads 0 \
  --vcf "${VCF}" \
  --bam "${BAM}" \
  --paired_end 1 \
  --mapq 255 --baseq 10 \
  --sample NA06986 \
  --chr 1 \
  --prefilter_hets 1 --reintegrate_vcf 1 \
  --o "${PHASER_PREFIX}" 2> "${PHASER_PREFIX}.time.txt"

test -s "${PHASER_PREFIX}.haplotypic_counts.txt"

echo "== phaser_gene_ae smoke test =="
FEATURES="${OUT_DIR}/features_chr1.bed"
printf "1\t0\t250000000\ttest_feature\n" > "${FEATURES}"

GENE_AE_DIR="${OUT_DIR}/gene_ae"
mkdir -p "${GENE_AE_DIR}"
GENE_AE_OUT="${GENE_AE_DIR}/NA06986.chr1.gene_ae.txt"

python3 "${ROOT}/phaser_gene_ae/phaser_gene_ae.py" \
  --haplotypic_counts "${PHASER_PREFIX}.haplotypic_counts.txt" \
  --features "${FEATURES}" \
  --o "${GENE_AE_OUT}"

test -s "${GENE_AE_OUT}"

echo "== phaser_pop smoke test (phaser_expr_matrix.py) =="
(
  cd "${OUT_DIR}"
  python3 "${ROOT}/phaser_pop/phaser_expr_matrix.py" \
    --gene_ae_dir "${GENE_AE_DIR}" \
    --features "${FEATURES}" \
    --t 1 \
    --o "${OUT_DIR}/expr_matrix"
)

test -s "${OUT_DIR}/expr_matrix.bed.gz"
test -s "${OUT_DIR}/expr_matrix.gw_phased.bed.gz"

echo "Smoke tests passed. Outputs are in: ${OUT_DIR}"


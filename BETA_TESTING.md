# phASER Fast Beta - Tester Guide

This bundle contains an optimized phASER build focused on speeding up the default SNP-only workflow (`--include_indels 0`).

## Quick Start (Recommended: conda)

1) Create the environment:

```bash
conda env create -f environment.yml
conda activate phaser-fast-beta
```

2) Download the official NA06986 test case (large: BAM is multiple GB):

```bash
./scripts/download_test_case.sh
```

3) Run the test case (single thread, default MAPQ=255):

```bash
./scripts/run_test_case.sh
```

Outputs go to `out/`.

## Alternative Setup (venv + pip)

If you already have external tools installed (`samtools`, `bcftools`, `bgzip`, `tabix`, `bedtools`) and working on your `PATH`:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Running On Your Own Data

From the bundle root:

```bash
python3 phaser/phaser.py \
  --threads 1 --io_threads 0 \
  --vcf /path/to/sample.vcf.gz \
  --bam /path/to/sample.bam \
  --paired_end 1 \
  --mapq 255 --baseq 10 \
  --sample SAMPLE_ID_IN_VCF \
  --o out/sample
```

Notes:
- `--mapq 255` is the default used in the upstream NA06986 example.
- `--include_indels 1` is supported but not part of the optimized path in this beta.
- This beta build is safe to run from paths containing spaces/parentheses (e.g. `Dropbox (Personal)`).

## Smoke Testing Other Tools

This bundle includes:
- `phaser_gene_ae/` (gene-level haplotypic expression)
- `phaser_pop/` (phASER-pop utilities)

You can run:

```bash
./scripts/smoke_test.sh
```

## Reporting Issues

Please include:
- Command used
- `*.log` output (phASER writes `<prefix>.log`)
- Platform (`uname -a`) and python version (`python3 --version`)
- Whether you used conda or pip

Known non-issue:
- `*.haplotypic_counts.txt` may show a different ordering of read-index lists (`aReads`/`bReads`) while counts remain the same.

# phASER Fast Beta (2026-02-03) - Executive Summary

This beta build targets phASER runtime. The goal was meaningful single-thread speedups (>=2x) without changing results for the default mode (`--include_indels 0`).

## What Changed (High Level)

1. **Eliminated the SAM-text bottleneck in read-to-variant mapping (biggest win).**
   - Old pipeline: `samtools view ... -> SAM text -> python parses SAM + CIGAR -> allele calls`
   - New pipeline (default): `samtools view -u ... -> BAM stream -> python reads BAM via pysam -> allele calls`

2. **Added a new SNP-only fast mapper for the default configuration (`--include_indels 0`).**
   - New script: `phaser/call_read_variant_map_bam_fast.py`
   - Uses `pysam` + CIGAR-aware reference mapping to extract alleles efficiently.
   - **Compatibility note:** preserves legacy behavior for spliced reads by *not* appending insertions that occur after a ref-skip (`N`) in the CIGAR (required for output parity).

3. **Kept an indel-capable fallback mapper.**
   - `--include_indels 1` routes to `phaser/call_read_variant_map_bam.py` (BAM-stream version of the legacy mapper logic).
   - Indels are still slower and are not the focus of this beta build.

4. **Made VCF output much faster (streaming + less repeated parsing).**
   - Streams input VCF via `tabix`/`gunzip` and writes output directly into `bgzip` (no large intermediate uncompressed `.vcf`).
   - Caches FORMAT parsing and avoids repeated expensive operations.

5. **Improved robustness for common macOS paths.**
   - Quoted shell-invoked file paths so running from directories with spaces/parentheses works (e.g. `Dropbox (Personal)`).
   - `phaser_pop/phaser_expr_matrix.py` and `phaser_pop/phaser_cis_var.py` were also updated for shell-safe paths; `phaser_cis_var.py` was made Python 3 friendly (StringIO + bytes/str fixes).

## Performance (Single Thread)

**Benchmark dataset:** NA06986 test case referenced in `phaser/README.md` (hg19).  
**Benchmark command (single thread, default MAPQ=255):**

```bash
python3 phaser/phaser.py \
  --threads 1 --io_threads 0 \
  --vcf data/NA06986.vcf.gz \
  --bam data/NA06986.2.M_111215_4.bam \
  --paired_end 1 --mapq 255 --baseq 10 \
  --sample NA06986 \
  --o out/phaser_test_case
```

**Measured on macOS 15.7.3 (arm64):**
- Baseline (original): ~635s
- Fast beta: ~76s
- **Speedup:** ~8.3x

## Output Compatibility

Compared to the original version on the NA06986 test case, outputs match (order-insensitive) for:
- `*.allele_config.txt`
- `*.allelic_counts.txt`
- `*.variant_connections.txt`
- `*.haplotypes.txt`

`*.haplotypic_counts.txt` matches on all count/stat columns; the read-index lists (`aReads`/`bReads`) may differ in ordering only.

## Known Limitations / Notes

- **Indels:** `--include_indels 1` uses the slower fallback mapper and is not optimized in this beta.
- **External tools still required:** `samtools`, `bcftools`, `bgzip`, `tabix`, `bedtools`.
- **Python 2 legacy:** `phaser_annotate` is legacy (python2-era) and was not modernized as part of this speedup work.

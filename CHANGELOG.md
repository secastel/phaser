# Changelog (Beta)

## 2026-02-03 - phaser_fast_beta_20260203

- Default read-to-variant mapping now uses a fast BAM-stream SNP-only mapper (`phaser/call_read_variant_map_bam_fast.py`) when `--include_indels 0`.
- `--include_indels 1` uses a BAM-stream version of the legacy mapper (`phaser/call_read_variant_map_bam.py`).
- VCF writing path is now streaming (no intermediate uncompressed VCF), with reduced per-record parsing overhead.
- Shell-invoked paths are quoted to support directories with spaces/parentheses (common on macOS).
- `phaser_pop/phaser_expr_matrix.py` and `phaser_pop/phaser_cis_var.py` updated for Python 3 + shell-safe paths.

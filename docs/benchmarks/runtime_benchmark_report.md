# Runtime Benchmark Report

This release includes a runtime comparison plot for baseline phASER vs tuned phASER modes.

## Datasets
- **NA06986** (as included in package manual)
- **NA20808** (as downloaded from 1KG + SHAPEIT5)

## Runtime Summary (wall-clock, 1 thread)

| Dataset | Baseline | Tuned (pref OFF) | Tuned (pref ON, reint ON) | Tuned (pref ON, reint OFF) |
|---|---:|---:|---:|---:|
| NA06986 | 641.62s | 88.01s | 89.05s | 87.37s |
| NA20808 | 1133.73s | 552.45s | 466.30s | 112.71s |

## Speedup vs Baseline

| Dataset | pref OFF | pref ON, reint ON | pref ON, reint OFF |
|---|---:|---:|---:|
| NA06986 | 7.29x | 7.21x | 7.34x |
| NA20808 | 2.05x | 2.43x | 10.06x |

## Notes
- `pref ON + reint OFF` is fastest when downstream only needs phased het-site output.
- `pref ON + reint ON` preserves backward-compatible full-VCF output.
- For cohorts where you must retain hom sites, run phASER with reintegration OFF and perform reintegration once as a separate batch step outside phASER to reduce total runtime.
- If `prefilter_hets` is off, reintegration is automatically disabled.

---
Prepared by [PEJ Lab](https://pejlab.org).

## Test environments
* local macOS Sequoia 15.0.1, R 4.6.0
* win-builder (devel and release)
* R-hub (Ubuntu, Windows, macOS)

## R CMD check results
0 errors | 0 warnings | 0 notes

## Resubmission notes
This is version 0.4.0. Major changes since 0.3.x:
* Added permutation tests for global and local SSIM in `ssim_polygon()` and `ssim_raster()`, with optional Benjamini–Hochberg FDR correction for local p-values.
* Added a new `transform` argument to `ssim_polygon()` and `ssim_bandwidth()` supporting `"normal_score"` (Blom normal scores), `"percentile"` (empirical percentiles), `"minmax"` (min–max normalization), and `"none"`.
* Redesigned `ssim_bandwidth()` with bias-variance trade-off plots and a new `option` argument (`"midpoint"`, `"lower"`, `"upper"`) for selecting bandwidth from the suggested range.
* Added new example data: BC Fire Weather Index rasters (`fwi_0816_bc`, `fwi_0818_bc`, `fwi_1101_bc`) and `Toronto_SSIM`.
* Removed `ssim_constant()` (constants `k1` and `k2` are now arguments of `ssim_polygon()`).

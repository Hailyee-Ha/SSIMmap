# SSIMmap 0.4.0

## Major changes
* `ssim_bandwidth()`: redesigned with `transform`, `option` arguments and bias-variance trade-off plots.
* `ssim_polygon()` and `ssim_raster()`: added permutation testing (`do_test`, `R`, `fdr`, `alpha`, `seed`) and FDR-corrected local p-values.

## Breaking changes
* Removed `ssim_constant()`. Constants `k1` and `k2` are now arguments of `ssim_polygon()`.

## New data
* Added `Toronto_SSIM`, `fwi_0816_bc`, `fwi_0818_bc`, `fwi_1101_bc`.

## New data
* Added `Toronto_SSIM` (sf polygon object).
* Added `fwi_0816_bc`, `fwi_0818_bc`, `fwi_1101_bc` (packed terra rasters of BC Fire Weather Index for Aug 16, Aug 18, and Nov 1, 2023).

#' Toronto neighbourhood deprivation indices and commuting pattern (2016)
#'
#' An \code{sf} polygon object for the Toronto Census Subdivision, Ontario,
#' Canada. Each polygon corresponds to a Census Tract (CT) and includes two
#' neighbourhood-level deprivation indices and one census variable, which can
#' be used as inputs to \code{\link{ssim_polygon}} and
#' \code{\link{ssim_bandwidth}}.
#'
#' @format An \code{sf} data frame with one row per Census Tract and the
#'   following key columns:
#'   \describe{
#'     \item{CIMD}{Canadian Index of Multiple Deprivation (situational
#'       dimension), 2016. Higher values indicate greater deprivation.}
#'     \item{PP}{Pampalon deprivation index (social dimension), 2016.}
#'     \item{Commute}{Percentage of households commuting within the
#'       Census Subdivision of residence, 2016.}
#'     \item{geometry}{\code{sfc} polygon geometry column.}
#'   }
#'
#' @source
#'   Statistics Canada, 2016 Census of Population.
#'
#' @examples
#' data("Toronto_SSIM")
#' plot(Toronto_SSIM["CIMD"], border = NA)
"Toronto_SSIM"


#' British Columbia Fire Weather Index — 16 August 2023
#'
#' Daily Canadian Fire Weather Index (FWI) for the province of British
#' Columbia, Canada, on 16 August 2023, at a 2 km spatial resolution. This
#' date represents conditions during the peak summer fire season, when fire
#' danger levels in British Columbia are typically high.
#'
#' The object is stored as a packed \code{terra} raster (created with
#' \code{terra::wrap()}) so that it can be saved in an R \code{.rda} file.
#' Use \code{terra::unwrap()} to convert it back to a regular \code{terra}
#' raster before passing it to \code{\link{ssim_raster}} or to other
#' \code{terra} functions.
#'
#' @format A packed \code{terra} raster object with one numeric layer
#'   (FWI values).
#'
#' @source
#'   Canadian Forest Fire Weather Index System; FWI values are derived from
#'   meteorological inputs (temperature, wind speed, relative humidity,
#'   precipitation).
#'
#' @examples
#' data("fwi_0816_bc")
#' fwi_0816_bc <- terra::unwrap(fwi_0816_bc)
#' terra::plot(fwi_0816_bc, main = "FWI - 16 August 2023")
"fwi_0816_bc"


#' British Columbia Fire Weather Index — 18 August 2023
#'
#' Daily Canadian Fire Weather Index (FWI) for the province of British
#' Columbia, Canada, on 18 August 2023, at a 2 km spatial resolution. This
#' date also represents conditions during the peak summer fire season and is
#' typically expected to show high spatial similarity to
#' \code{\link{fwi_0816_bc}}.
#'
#' The object is stored as a packed \code{terra} raster (created with
#' \code{terra::wrap()}). Use \code{terra::unwrap()} to convert it back to a
#' regular \code{terra} raster before use.
#'
#' @format A packed \code{terra} raster object with one numeric layer
#'   (FWI values).
#'
#' @source Canadian Forest Fire Weather Index System.
#'
#' @examples
#' data("fwi_0818_bc")
#' fwi_0818_bc <- terra::unwrap(fwi_0818_bc)
#' terra::plot(fwi_0818_bc, main = "FWI - 18 August 2023")
"fwi_0818_bc"


#' British Columbia Fire Weather Index — 1 November 2023
#'
#' Daily Canadian Fire Weather Index (FWI) for the province of British
#' Columbia, Canada, on 1 November 2023, at a 2 km spatial resolution. This
#' date represents late-fall conditions, when fire danger across most of
#' British Columbia is typically much lower than during the summer season due
#' to cooler temperatures and increased precipitation. The spatial structure
#' of this map is therefore expected to differ substantially from the summer
#' FWI maps (\code{\link{fwi_0816_bc}}, \code{\link{fwi_0818_bc}}).
#'
#' The object is stored as a packed \code{terra} raster (created with
#' \code{terra::wrap()}). Use \code{terra::unwrap()} to convert it back to a
#' regular \code{terra} raster before use.
#'
#' @format A packed \code{terra} raster object with one numeric layer
#'   (FWI values).
#'
#' @source Canadian Forest Fire Weather Index System.
#'
#' @examples
#' data("fwi_1101_bc")
#' fwi_1101_bc <- terra::unwrap(fwi_1101_bc)
#' terra::plot(fwi_1101_bc, main = "FWI - 1 November 2023")
"fwi_1101_bc"

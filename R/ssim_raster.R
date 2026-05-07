#' SSIM index for raster images
#'
#' This function calculates the Structural Similarity Index Measure (SSIM)
#' between two raster images. It can return either global SSIM summaries or
#' local SSIM raster layers. Optional permutation-based statistical inference
#' can also be performed for global or local SSIM values.
#'
#' @param map1 A single-layer \code{terra::SpatRaster} object representing the first raster map.
#' @param map2 A single-layer \code{terra::SpatRaster} object representing the second raster map.
#' @param global Logical. If \code{TRUE}, the function returns global summary
#'   statistics of SSIM, SIM, SIV, and SIP. If \code{FALSE}, it returns local
#'   raster layers for each cell. Default is \code{TRUE}.
#' @param w Integer. Radius of the local moving window used to calculate local
#'   SSIM components. For example, \code{w = 1} uses a 3 x 3 window, and
#'   \code{w = 2} uses a 5 x 5 window. Default is \code{1}.
#' @param transform Character. Transformation applied to both raster maps before
#'   SSIM calculation. Options are \code{"normal_score"}, \code{"percentile"},
#'   \code{"none"}, and \code{"minmax"}. Default is \code{"normal_score"}.
#' @param k1 Numeric. Constant used to stabilize the luminance component of SSIM.
#'   If \code{NULL}, the default value \code{0.01} is used.
#' @param k2 Numeric. Constant used to stabilize the contrast and structure
#'   components of SSIM. If \code{NULL}, the default value \code{0.03} is used.
#' @param do_test Logical. If \code{TRUE}, a permutation-based statistical test
#'   is performed. Default is \code{FALSE}.
#' @param local_test Logical. If \code{TRUE} and \code{global = FALSE}, local
#'   permutation tests are performed for each raster cell. If \code{FALSE},
#'   local SSIM layers are returned without p-values. Default is \code{FALSE}.
#' @param R Integer. Number of permutations used for statistical inference.
#'   Default is \code{1000}.
#' @param fdr Logical. If \code{TRUE}, p-values are adjusted using the
#'   Benjamini-Hochberg false discovery rate correction. Default is \code{TRUE}.
#' @param alpha Numeric. Significance level used to determine significant SSIM
#'   values. Default is \code{0.05}.
#' @param seed Optional integer. Random seed used for permutation testing.
#'   Default is \code{NULL}.
#'
#' @return
#' If \code{global = TRUE} and \code{do_test = FALSE}, a list containing a summary
#' table of global SSIM, SIM, SIV, and SIP values is returned.
#'
#' If \code{global = TRUE} and \code{do_test = TRUE}, a list containing the
#' summary table, permutation null distributions, p-values, and optionally
#' q-values is returned.
#'
#' If \code{global = FALSE} and \code{do_test = FALSE}, a
#' \code{terra::SpatRaster} object with four layers is returned:
#' \code{SSIM}, \code{SIM}, \code{SIV}, and \code{SIP}.
#'
#' If \code{global = FALSE}, \code{do_test = TRUE}, and
#' \code{local_test = TRUE}, a \code{terra::SpatRaster} object is returned with
#' local SSIM component layers, p-value layers, optionally q-value layers, and
#' significance layers.
#'
#' @details
#' The SSIM index is calculated as the product of three components:
#' similarity in mean intensity, similarity in variance, and similarity in
#' spatial pattern. These components are returned as \code{SIM}, \code{SIV}, and
#' \code{SIP}, respectively.
#'
#' The local SSIM calculation uses a square moving window with size
#' \code{(2 * w + 1) x (2 * w + 1)}.
#'
#' When \code{do_test = TRUE}, permutation testing is performed by randomly
#' permuting the values of \code{map2} over the overlapping non-\code{NA} cells
#' and recalculating SSIM values.
#'
#' @examples
#' \dontrun{
#' library(terra)
#'
#' # Global SSIM summary
#' g1 <- ssim_raster(
#'   map1,
#'   map2,
#'   global = TRUE,
#'   w = 1
#' )
#'
#' # Local SSIM raster layers
#' l1 <- ssim_raster(
#'   map1,
#'   map2,
#'   global = FALSE,
#'   w = 1
#' )
#'
#' # Local SSIM with permutation test
#' l2 <- ssim_raster(
#'   map1,
#'   map2,
#'   global = FALSE,
#'   w = 1,
#'   do_test = TRUE,
#'   local_test = TRUE,
#'   R = 999
#' )
#' }
#'
#' @importFrom terra global focal crop mask compareGeom values nlyr ifel
#' @export
ssim_raster <- function(map1,
                        map2,
                        global = TRUE,
                        w = 1,
                        transform = c("normal_score", "percentile", "none", "minmax"),
                        k1 = NULL,
                        k2 = NULL,
                        do_test = FALSE,
                        local_test = FALSE,
                        R = 1000,
                        fdr = TRUE,
                        alpha = 0.05,
                        seed = NULL) {

  transform <- match.arg(transform)

  if (!inherits(map1, "SpatRaster") || !inherits(map2, "SpatRaster")) {
    stop("map1 and map2 must be terra SpatRaster objects.")
  }

  if (terra::nlyr(map1) != 1 || terra::nlyr(map2) != 1) {
    stop("map1 and map2 must each have exactly one layer.")
  }

  if (!terra::compareGeom(map1, map2, stopOnError = FALSE)) {
    stop("map1 and map2 must have the same extent, resolution, CRS, and geometry.")
  }

  if (!is.logical(global) || length(global) != 1) {
    stop("global must be TRUE or FALSE.")
  }

  if (!is.logical(do_test) || length(do_test) != 1) {
    stop("do_test must be TRUE or FALSE.")
  }

  if (!is.logical(local_test) || length(local_test) != 1) {
    stop("local_test must be TRUE or FALSE.")
  }

  if (w < 1 || w != floor(w)) {
    stop("w must be a positive integer. w = 1 means a 3x3 window.")
  }

  if (do_test && (R < 1 || R != floor(R))) {
    stop("R must be a positive integer.")
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  map1 <- map1[[1]]
  map2 <- map2[[1]]

  # -------------------------------------------------------------------
  # Helper functions
  # -------------------------------------------------------------------

  .minmax <- function(x) {
    mn <- terra::global(x, "min", na.rm = TRUE)[1, 1]
    mx <- terra::global(x, "max", na.rm = TRUE)[1, 1]
    c(min = mn, max = mx)
  }

  .rank_transform <- function(x, type = c("percentile", "normal_score")) {
    type <- match.arg(type)

    vals <- terra::values(x, mat = FALSE)
    idx <- which(!is.na(vals))

    out <- rep(NA_real_, length(vals))

    if (length(idx) == 0) {
      y <- x
      terra::values(y) <- out
      return(y)
    }

    rnk <- rank(vals[idx], ties.method = "average")
    p <- (rnk - 0.5) / length(rnk)

    if (type == "percentile") {
      out[idx] <- p
    }

    if (type == "normal_score") {
      out[idx] <- stats::qnorm(p)
    }

    y <- x
    terra::values(y) <- out
    y
  }

  .scale01 <- function(x) {
    mm <- .minmax(x)

    if (!is.finite(mm["min"]) || !is.finite(mm["max"])) {
      stop("Raster has no finite values.")
    }

    if (mm["max"] == mm["min"]) {
      return(x * 0)
    }

    (x - mm["min"]) / (mm["max"] - mm["min"])
  }

  .apply_transform <- function(x, transform) {
    if (transform == "none") {
      return(x)
    }

    if (transform == "minmax") {
      return(.scale01(x))
    }

    if (transform == "percentile") {
      return(.rank_transform(x, "percentile"))
    }

    if (transform == "normal_score") {
      return(.rank_transform(x, "normal_score"))
    }
  }

  .ssim_layers <- function(x, y, w, k1, k2) {
    if (is.null(k1)) k1 <- 0.01
    if (is.null(k2)) k2 <- 0.03

    if (k1 <= 0 || k1 >= 1) {
      stop("k1 must be greater than 0 and smaller than 1.")
    }

    if (k2 <= 0 || k2 >= 1) {
      stop("k2 must be greater than 0 and smaller than 1.")
    }

    mm_x <- .minmax(x)
    mm_y <- .minmax(y)

    data_min <- min(mm_x["min"], mm_y["min"], na.rm = TRUE)
    data_max <- max(mm_x["max"], mm_y["max"], na.rm = TRUE)
    L_range <- data_max - data_min

    if (!is.finite(L_range) || L_range == 0) {
      L_range <- 1
    }

    C1 <- (k1 * L_range)^2
    C2 <- (k2 * L_range)^2
    C3 <- C2 / 2

    win_size <- (2 * w) + 1

    kernel <- matrix(
      1 / (win_size^2),
      nrow = win_size,
      ncol = win_size
    )

    fmean <- function(z) {
      terra::focal(
        z,
        w = kernel,
        fun = "sum",
        na.policy = "all",
        fillvalue = NA
      )
    }

    mu1 <- fmean(x)
    mu2 <- fmean(y)

    sig1 <- sqrt(abs(fmean(x * x) - mu1 * mu1))
    sig2 <- sqrt(abs(fmean(y * y) - mu2 * mu2))
    sig12 <- fmean(x * y) - mu1 * mu2

    SIM <- ((2 * mu1 * mu2) + C1) / ((mu1^2) + (mu2^2) + C1)
    SIV <- ((2 * sig1 * sig2) + C2) / ((sig1^2) + (sig2^2) + C2)
    SIP <- (sig12 + C3) / ((sig1 * sig2) + C3)

    SSIM <- SIM * SIV * SIP

    out <- c(SSIM, SIM, SIV, SIP)
    names(out) <- c("SSIM", "SIM", "SIV", "SIP")
    out
  }

  .summary_table <- function(x) {
    data.frame(
      metric = names(x),
      mean = as.numeric(terra::global(x, "mean", na.rm = TRUE)[, 1]),
      min  = as.numeric(terra::global(x, "min",  na.rm = TRUE)[, 1]),
      max  = as.numeric(terra::global(x, "max",  na.rm = TRUE)[, 1]),
      sd   = as.numeric(terra::global(x, stats::sd, na.rm = TRUE)[, 1]),
      row.names = NULL
    )
  }

  .permute_raster <- function(x, valid_idx) {
    vals <- terra::values(x, mat = FALSE)

    out <- rep(NA_real_, length(vals))
    out[valid_idx] <- sample(vals[valid_idx], length(valid_idx), replace = FALSE)

    y <- x
    terra::values(y) <- out
    y
  }

  .p_adjust_raster <- function(p_raster) {
    q_raster <- p_raster

    for (i in seq_len(terra::nlyr(p_raster))) {
      vals <- terra::values(p_raster[[i]], mat = FALSE)
      idx <- which(!is.na(vals))

      q_vals <- rep(NA_real_, length(vals))
      q_vals[idx] <- stats::p.adjust(vals[idx], method = "BH")

      terra::values(q_raster[[i]]) <- q_vals
    }

    names(q_raster) <- sub("_p$", "_q", names(p_raster))
    q_raster
  }

  .sig_raster <- function(test_raster, alpha) {
    sig <- test_raster <= alpha
    sig <- terra::ifel(is.na(sig), NA, sig)
    names(sig) <- sub("_(p|q)$", "_sig", names(test_raster))
    sig
  }

  # -------------------------------------------------------------------
  # Pre-processing
  # -------------------------------------------------------------------

  valid <- !is.na(map1) & !is.na(map2)

  map1 <- terra::mask(map1, valid, maskvalues = FALSE)
  map2 <- terra::mask(map2, valid, maskvalues = FALSE)

  map1 <- .apply_transform(map1, transform)
  map2 <- .apply_transform(map2, transform)

  valid_vals <- terra::values(valid, mat = FALSE)
  valid_idx <- which(!is.na(valid_vals) & valid_vals == 1)

  if (length(valid_idx) == 0) {
    stop("There are no overlapping non-NA cells between map1 and map2.")
  }

  # -------------------------------------------------------------------
  # Observed SSIM
  # -------------------------------------------------------------------

  observed <- .ssim_layers(map1, map2, w, k1, k2)
  observed <- terra::mask(observed, valid, maskvalues = FALSE)

  obs_summary <- .summary_table(observed)

  # -------------------------------------------------------------------
  # Case 1: No statistical test
  # -------------------------------------------------------------------

  if (!do_test) {
    if (global) {
      return(list(
        summary = obs_summary,
        settings = list(
          global = global,
          w = w,
          transform = transform,
          do_test = do_test,
          local_test = local_test
        )
      ))
    }

    return(observed)
  }

  # -------------------------------------------------------------------
  # Case 2: global = TRUE, do_test = TRUE
  # Global permutation test
  # -------------------------------------------------------------------

  if (global) {
    if (local_test) {
      warning("local_test = TRUE is ignored because global = TRUE.")
    }

    null_means <- matrix(
      NA_real_,
      nrow = R,
      ncol = terra::nlyr(observed)
    )

    colnames(null_means) <- names(observed)

    for (i in seq_len(R)) {
      map2_perm <- .permute_raster(map2, valid_idx)
      perm_layers <- .ssim_layers(map1, map2_perm, w, k1, k2)
      perm_layers <- terra::mask(perm_layers, valid, maskvalues = FALSE)

      null_means[i, ] <- as.numeric(
        terra::global(perm_layers, "mean", na.rm = TRUE)[, 1]
      )
    }

    obs_means <- obs_summary$mean

    p_values <- vapply(seq_along(obs_means), function(j) {
      (sum(null_means[, j] >= obs_means[j], na.rm = TRUE) + 1) / (R + 1)
    }, numeric(1))

    obs_summary$p_value <- p_values

    if (isTRUE(fdr)) {
      obs_summary$q_value <- stats::p.adjust(obs_summary$p_value, method = "BH")
      obs_summary$significant <- obs_summary$q_value <= alpha
    } else {
      obs_summary$significant <- obs_summary$p_value <= alpha
    }

    return(list(
      summary = obs_summary,
      null_distribution = as.data.frame(null_means),
      settings = list(
        global = global,
        w = w,
        transform = transform,
        do_test = do_test,
        local_test = local_test,
        R = R,
        fdr = fdr,
        alpha = alpha,
        seed = seed
      )
    ))
  }

  # -------------------------------------------------------------------
  # Case 3: global = FALSE, do_test = TRUE, local_test = FALSE
  # Return local SSIM only, no local p-value
  # -------------------------------------------------------------------

  if (!global && do_test && !local_test) {
    warning(
      "do_test = TRUE but local_test = FALSE. ",
      "Returning local SSIM layers only. ",
      "Use local_test = TRUE for local p-value/q-value layers."
    )

    return(observed)
  }

  # -------------------------------------------------------------------
  # Case 4: global = FALSE, do_test = TRUE, local_test = TRUE
  # Local permutation test
  # -------------------------------------------------------------------

  count_ge <- observed * 0
  count_ge <- terra::ifel(is.na(count_ge), 0, count_ge)

  for (i in seq_len(R)) {
    map2_perm <- .permute_raster(map2, valid_idx)
    perm_layers <- .ssim_layers(map1, map2_perm, w, k1, k2)
    perm_layers <- terra::mask(perm_layers, valid, maskvalues = FALSE)

    ge <- perm_layers >= observed
    ge <- terra::ifel(is.na(ge), 0, ge)

    count_ge <- count_ge + ge
  }

  p_layers <- (count_ge + 1) / (R + 1)
  p_layers <- terra::mask(p_layers, observed[[1]])
  names(p_layers) <- paste0(names(observed), "_p")

  if (isTRUE(fdr)) {
    q_layers <- .p_adjust_raster(p_layers)
    sig_layers <- .sig_raster(q_layers, alpha)

    out <- c(observed, p_layers, q_layers, sig_layers)
  } else {
    sig_layers <- .sig_raster(p_layers, alpha)

    out <- c(observed, p_layers, sig_layers)
  }

  return(out)
}

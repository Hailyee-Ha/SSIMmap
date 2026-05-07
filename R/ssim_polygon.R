#' Structural similarity index (SSIM) for polygon maps
#'
#' Computes local and global SSIM (and its components SIM, SIV, SIP) for two
#' polygon maps using an adaptive k-nearest-neighbour (k-NN) Gaussian kernel.
#' Input variables can be optionally transformed (e.g. rank-based inverse
#' normal scores or min–max normalization). 
#' Optionally performs permutation tests for global and local significance with BH-FDR correction.
#' Under the null hypothesis for the permutation test, both variables are
#' randomly permuted at each iteration, breaking spatial structure and
#' association between them. Local p-values can be optionally adjusted using
#' the Benjamini–Hochberg FDR procedure.
#'
#' @param shape An \code{sf} polygon object.
#' @param map1 Character string; column name in \code{shape} for the first map.
#' @param map2 Character string; column name in \code{shape} for the second map.
#' @param global Logical. If \code{TRUE}, compute and print a summary table of
#'   global SSIM, SIM, SIV, SIP. If \code{FALSE}, return an \code{sf} object
#'   with local metrics (and local p-/q-values if \code{do_test = TRUE}).
#' @param bandwidth Integer; adaptive k-NN size (number of neighbours). The
#'   default is \code{ceiling(sqrt(n))}, where \code{n} is the number of
#'   polygons. Must be at least 3 and not exceed \code{n}.
#' @param transform One of \code{c("normal_score", "percentile", "none", "minmax")}.
#'   \code{"normal_score"} applies a Blom normal scores transform;
#'   \code{"percentile"} maps values to empirical percentiles in (0, 1);
#'   \code{"minmax"} applies min–max normalisation to [0, 1]; and
#'   \code{"none"} leaves the variables on their original scale.
#' @param k1,k2 SSIM constants. If \code{NULL}, defaults are \code{k1 = 0.01},
#'   \code{k2 = 0.03}.
#' @param do_test Logical; if \code{TRUE}, perform permutation tests for global
#'   and local SSIM. Local p-values can be FDR-adjusted if \code{fdr = TRUE}.
#' @param R Integer; number of permutations for the significance tests
#'   (default \code{R = 1000}).
#' @param fdr Logical; if \code{TRUE} (default), apply Benjamini–Hochberg
#'   FDR correction to local p-values.
#' @param alpha Numeric; significance threshold for local results (default
#'   \code{0.05}). Local features with \code{q_value < alpha} are flagged as
#'   significant.
#' @param seed Optional integer; random seed for reproducibility of the
#'   permutation tests.
#'
#' @return
#'   If \code{global = TRUE} and \code{do_test = FALSE}, the function prints a
#'   knitr table summarising the global mean, minimum, maximum, and standard
#'   deviation of SSIM, SIM, SIV, and SIP, and returns this \code{data.frame}
#'   (invisibly).
#'
#'   If \code{global = TRUE} and \code{do_test = TRUE}, the function prints the
#'   same summary table plus a table of global permutation p-values (two-sided)
#'   for SSIM, SIM, SIV, and SIP, and returns (invisibly) a \code{list} with
#'   components:
#'   \itemize{
#'     \item \code{summary}: \code{data.frame} with global SSIM/SIM/SIV/SIP
#'       summary statistics;
#'     \item \code{p_global}: \code{data.frame} with global means and
#'       permutation p-values for SSIM, SIM, SIV, SIP.
#'   }
#'
#'   If \code{global = FALSE}, the function returns an \code{sf} object equal to
#'   \code{shape} with additional columns:
#'   \itemize{
#'     \item \code{SSIM}, \code{SIM}, \code{SIV}, \code{SIP}: local similarity
#'       metrics for each polygon;
#'     \item \code{p_value}, \code{q_value}, \code{sig} (only if
#'       \code{do_test = TRUE}): local permutation p-value, FDR-adjusted q-value,
#'       and a logical flag indicating significance \code{(q_value < alpha)}.
#'   }
#'
#' @importFrom sf st_geometry st_centroid st_coordinates
#' @importFrom FNN get.knn
#' @importFrom stats qnorm sd quantile p.adjust
#' @importFrom knitr kable
#'
#' @export
#'
ssim_polygon <- function(
    shape, map1, map2,
    global = TRUE,
    bandwidth = NULL,
    transform = c("normal_score","percentile","none","minmax"),
    k1 = NULL, k2 = NULL,
    do_test = FALSE, R = 1000, fdr = TRUE, alpha = 0.05,
    seed = NULL
){
  transform <- match.arg(transform)
  if (identical(map1, map2)) stop("variables are identical")
  if (!map1 %in% names(shape) || !map2 %in% names(shape))
    stop("map1/map2 must be column names in 'shape'")
  
  n <- nrow(shape)
  if (is.null(bandwidth)) bandwidth <- ceiling(sqrt(n))
  if (bandwidth < 3) stop("bandwidth must be >= 3")
  if (bandwidth > n) stop("bandwidth cannot exceed number of polygons")
  
  if (!is.null(seed)) set.seed(seed)
  
  # --- helpers: transforms ---
  normal_scores <- function(x, ties="average"){
    # Blom normal scores transformation
    r <- rank(x, ties.method = ties); n <- length(x)
    u <- (r - 3/8) / (n + 1/4)
    qnorm(u)
  }
  
  to_percentile <- function(x, ties="average"){
    # Map to (0,1) empirical percentiles
    rank(x, ties.method = ties)/(length(x)+1)
  }
  
  to_minmax <- function(x){
    # Standard min-max normalization to [0,1]
    rng <- range(x, na.rm = TRUE)
    if (!all(is.finite(rng)) || rng[1] == rng[2]) {
      stop("Min-max normalization is not defined when all values are identical or non-finite.")
    }
    (x - rng[1]) / (rng[2] - rng[1])
  }
  
  # --- extract & transform variables ---
  A_raw <- as.numeric(shape[[map1]])
  B_raw <- as.numeric(shape[[map2]])
  
  A <- switch(
    transform,
    normal_score = normal_scores(A_raw),
    percentile   = to_percentile(A_raw),
    none         = A_raw,
    minmax       = to_minmax(A_raw)
  )
  
  B <- switch(
    transform,
    normal_score = normal_scores(B_raw),
    percentile   = to_percentile(B_raw),
    none         = B_raw,
    minmax       = to_minmax(B_raw)
  )
  
  # --- geometry & neighbours (adaptive k-NN) ---
  coords <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(shape)))
  knn <- FNN::get.knn(coords, k = bandwidth)
  idx  <- knn$nn.index
  dist <- knn$nn.dist
  hi   <- apply(dist, 1, max)               # distance to k-th neighbour
  W    <- exp(-0.5 * (dist / hi)^2)         # Gaussian kernel
  Wn   <- W / rowSums(W)                    # row-normalize weights
  
  # --- local stats given A,B ---
  lw_stats_two <- function(A, B){
    mA <- rowSums(Wn * A[idx]); mB <- rowSums(Wn * B[idx])
    # local variances (expectation of squared deviations with normalized W)
    vA <- rowSums(Wn * (A[idx] - mA)^2)
    vB <- rowSums(Wn * (B[idx] - mB)^2)
    cAB<- rowSums(Wn * (A[idx] - mA) * (B[idx] - mB))
    list(mA=mA, mB=mB, vA=vA, vB=vB, cAB=cAB)
  }
  st_obs <- lw_stats_two(A, B)
  
  # --- SSIM constants (robust) ---
  if (is.null(k1)) k1 <- 0.01
  if (is.null(k2)) k2 <- 0.03
  if (k1 <= 0 || k2 <= 0 || k1 >= 1 || k2 >= 1) stop("k1/k2 must be in (0,1)")
  
  both <- c(A, B)
  Rrob <- as.numeric(
    quantile(both, 0.99, na.rm=TRUE) - quantile(both, 0.01, na.rm=TRUE)
  )
  if (!is.finite(Rrob) || Rrob <= 0) Rrob <- diff(range(both, na.rm=TRUE))
  C1 <- (k1 * Rrob)^2; C2 <- (k2 * Rrob)^2; C3 <- C2/2
  
  # --- components (raw) ---
  SIM_raw <- (2*st_obs$mA*st_obs$mB + C1) /
    (st_obs$mA^2 + st_obs$mB^2 + C1)
  SIV_raw <- (2*sqrt(st_obs$vA)*sqrt(st_obs$vB) + C2) /
    (st_obs$vA + st_obs$vB + C2)
  SIP_raw <- (st_obs$cAB + C3) /
    (sqrt(st_obs$vA)*sqrt(st_obs$vB) + C3)
  
  tol <- 1e-8
  
  # Hard range check: if any metric is clearly outside its theoretical range,
  # stop and ask the user to try a different transform option.
  check_range <- function(x, lo, hi, name){
    if (any(x < lo - tol | x > hi + tol, na.rm = TRUE)) {
      stop(
        sprintf(
          "%s values are outside their theoretical range [%g, %g]. This may indicate numerical instability. Please try a different 'transform' option.",
          name, lo, hi
        )
      )
    }
  }
  
  check_range(SIM_raw, 0, 1,  "SIM")
  check_range(SIV_raw, 0, 1,  "SIV")
  check_range(SIP_raw, -1, 1, "SIP")
  
  SIM <- SIM_raw
  SIV <- SIV_raw
  SIP <- SIP_raw
  SSIM <- SIM * SIV * SIP
  
  if (global && !do_test){
    res <- data.frame(
      Statistic = c("Mean","Min","Max","SD"),
      SSIM = c(mean(SSIM), min(SSIM), max(SSIM), sd(SSIM)),
      SIM  = c(mean(SIM ), min(SIM ), max(SIM ), sd(SIM )),
      SIV  = c(mean(SIV ), min(SIV ), max(SIV ), sd(SIV )),
      SIP  = c(mean(SIP ), min(SIP ), max(SIP ), sd(SIP ))
    )
    print(knitr::kable(res, row.names = FALSE))
    invisible(res)
  } else {
    p_loc <- q_loc <- sig <- NULL
    
    if (do_test){
      S_g_null   <- numeric(R)
      SIM_g_null <- numeric(R)
      SIV_g_null <- numeric(R)
      SIP_g_null <- numeric(R)
      S_i_null   <- matrix(NA_real_, nrow = n, ncol = R)
      
      for (r in seq_len(R)){
        
        # ---- PERMUTATION STEP (UPDATED) ----
        # Under the null hypothesis, both A and B are randomly permuted.
        # This breaks any spatial structure and association between them.
        Ap <- sample(A, replace = FALSE)
        Bp <- sample(B, replace = FALSE)
        
        st  <- lw_stats_two(Ap, Bp)
        
        SIMp_raw <- (2*st$mA*st$mB + C1) /
          (st$mA^2 + st$mB^2 + C1)
        SIVp_raw <- (2*sqrt(st$vA)*sqrt(st$vB) + C2) /
          (st$vA + st$vB + C2)
        SIPp_raw <- (st$cAB + C3) /
          (sqrt(st$vA)*sqrt(st$vB) + C3)
        
        # Range check also applied to permutation metrics
        check_range(SIMp_raw, 0, 1,  "SIM (permutation)")
        check_range(SIVp_raw, 0, 1,  "SIV (permutation)")
        check_range(SIPp_raw, -1, 1, "SIP (permutation)")
        
        SIMp <- SIMp_raw
        SIVp <- SIVp_raw
        SIPp <- SIPp_raw
        SSi  <- SIMp * SIVp * SIPp
        
        SIM_g_null[r] <- mean(SIMp)
        SIV_g_null[r] <- mean(SIVp)
        SIP_g_null[r] <- mean(SIPp)
        S_g_null[r]   <- mean(SSi)
        
        S_i_null[, r] <- SSi
      }
      
      # Observed global means
      SIM_g  <- mean(SIM); SIV_g <- mean(SIV); SIP_g <- mean(SIP); SSIM_g <- mean(SSIM)
      
      pval_two_sided <- function(null_vec, obs_mean){
        mu <- mean(null_vec)
        (sum(abs(null_vec - mu) >= abs(obs_mean - mu)) + 1) / (length(null_vec) + 1)
      }
      p_global_SSIM <- pval_two_sided(S_g_null,   SSIM_g)
      p_global_SIM  <- pval_two_sided(SIM_g_null, SIM_g)
      p_global_SIV  <- pval_two_sided(SIV_g_null, SIV_g)
      p_global_SIP  <- pval_two_sided(SIP_g_null, SIP_g)
      
      # Local SSIM p / q
      mu_i  <- rowMeans(S_i_null)
      p_loc <- (rowSums(abs(S_i_null - mu_i) >= abs(SSIM - mu_i)) + 1) / (R + 1)
      q_loc <- if (fdr) p.adjust(p_loc, method = "BH") else p_loc
      sig   <- q_loc < alpha
    }
    
    # ------ return branch ------
    if (global) {
      res <- data.frame(
        Statistic = c("Mean","Min","Max","SD"),
        SSIM = c(mean(SSIM), min(SSIM), max(SSIM), sd(SSIM)),
        SIM  = c(mean(SIM ), min(SIM ), max(SIM ), sd(SIM )),
        SIV  = c(mean(SIV ), min(SIV ), max(SIV ), sd(SIV )),
        SIP  = c(mean(SIP ), min(SIP ), max(SIP ), sd(SIP ))
      )
      print(knitr::kable(res, row.names = FALSE))
      
      if (do_test){
        res_p <- data.frame(
          Component   = c("SSIM","SIM","SIV","SIP"),
          Global_Mean = c(SSIM_g, SIM_g, SIV_g, SIP_g),
          P_value     = c(p_global_SSIM, p_global_SIM, p_global_SIV, p_global_SIP)
        )
        cat("\nGlobal permutation p-values (two-sided):\n")
        print(knitr::kable(res_p, row.names = FALSE, digits = 4))
        return(invisible(list(summary = res, p_global = res_p)))
      } else {
        return(invisible(res))
      }
      
    } else {
      # global == FALSE → always return sf with local metrics
      out_sf <- shape
      out_sf$SSIM <- SSIM
      out_sf$SIM  <- SIM
      out_sf$SIV  <- SIV
      out_sf$SIP  <- SIP
      if (do_test){
        out_sf$p_value <- p_loc
        out_sf$q_value <- q_loc
        out_sf$sig     <- sig
      }
      return(out_sf)
    }
  }
}

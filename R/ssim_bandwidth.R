#' Bandwidth selection for polygon maps
#'
#' Computes bias–variance trade-off curves for two polygon maps using an
#' adaptive Gaussian kernel. The two input variables can be optionally
#' transformed (e.g. normal scores or min–max normalisation) before smoothing.
#' A suggested bandwidth range is obtained by balancing standardized bias and
#' variance for each variable, and a single bandwidth is chosen from this range.
#'
#' This function does *not* compute SSIM; it is intended to be used in
#' combination with \code{ssim_polygon()}, which can then be called with the
#' selected bandwidth.
#'
#' @param shape An \code{sf} polygon object with attributes containing the two
#'   variables of interest.
#' @param map1 Character string; column name in \code{shape} for the first map.
#' @param map2 Character string; column name in \code{shape} for the second map.
#' @param max_bandwidth Integer (at least 12). Upper bound of the bandwidth
#'   \code{k} for the adaptive k-nearest-neighbour kernel. Must not exceed the
#'   number of polygons.
#' @param transform One of \code{c("normal_score","percentile","none","minmax")}.
#'   \code{"normal_score"} applies a Blom normal scores transform;
#'   \code{"percentile"} maps values to empirical percentiles in (0, 1);
#'   \code{"minmax"} applies min–max normalisation to [0, 1]; and
#'   \code{"none"} leaves the variables on their original scale.
#' @param option Character string specifying how to choose a single bandwidth
#'   from the bias–variance trade-off range: one of
#'   \code{"midpoint"}, \code{"lower"}, or \code{"upper"} (default
#'   \code{"midpoint"}). The square-root of the number of polygons (\code{sqrt(n)})
#'   is also shown on the plot as a reference.
#'
#' @return A list with components:
#'   \item{plot}{A \code{ggplot} object showing the bias–variance trade-off
#'     curves for both variables, the suggested bandwidth range, the chosen
#'     bandwidth, and \code{sqrt(n)} as a reference.}
#'   \item{bandwidth}{The chosen bandwidth \code{k_star}.}
#'   \item{tradeoff}{A list with elements
#'     \code{map1}, \code{map2}, \code{kA}, \code{kB}, and \code{k_sqrt},
#'     where \code{map1} and \code{map2} are data frames containing the
#'     standardized bias and variance across \code{k} for each variable.}
#'
#' @details
#' The kernel is Gaussian with adaptive bandwidth:
#' \deqn{w_{ij} = \exp\left\{-\tfrac12 (d_{ij} / h_i)^2 \right\},}
#' where \eqn{d_{ij}} is the distance between polygon centroids and
#' \eqn{h_i} is the distance from polygon \eqn{i} to its \eqn{k}-th nearest
#' neighbour. Bandwidth selection is based purely on bias–variance trade-off
#' and does not involve any SSIM computation or permutation testing.
#'
#' @importFrom sf st_geometry st_centroid st_coordinates
#' @importFrom FNN get.knn
#' @importFrom ggplot2 ggplot geom_line geom_rect geom_vline labs geom_text aes
#' @importFrom stats qnorm var
#' @importFrom scales rescale
#' @importFrom utils globalVariables
#' @export
ssim_bandwidth <- function(
    shape, map1, map2,
    max_bandwidth,
    transform = c("normal_score","percentile","none","minmax"),
    option = "midpoint"
){
  ## ---------- checks ----------
  transform <- match.arg(transform)
  if (map1 == map2) stop("variables are identical")
  if (max_bandwidth < 12) stop("max_bandwidth must be >= 12")
  n <- nrow(shape)
  if (n < max_bandwidth) stop("max_bandwidth is larger than number of polygons (n)")

  ## ---------- helpers ----------
  normal_scores <- function(x, ties = "average"){
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

  ## ---------- extract & transform ----------
  dat <- as.data.frame(shape)
  xA <- dat[[map1]]
  xB <- dat[[map2]]

  xA_t <- switch(
    transform,
    normal_score = normal_scores(xA),
    percentile   = to_percentile(xA),
    none         = xA,
    minmax       = to_minmax(xA)
  )

  xB_t <- switch(
    transform,
    normal_score = normal_scores(xB),
    percentile   = to_percentile(xB),
    none         = xB,
    minmax       = to_minmax(xB)
  )

  ## ---------- geometry & k-NN (for all k we reuse the same neighbour order) ----------
  coords <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(shape)))
  nn <- FNN::get.knn(coords, k = max_bandwidth)

  ## ---------- local weighted mean for a given k ----------
  lw_mean <- function(x, k){
    idx  <- nn$nn.index[, 1:k, drop=FALSE]
    dist <- nn$nn.dist[,  1:k, drop=FALSE]
    hi   <- apply(dist, 1, max)
    W    <- exp(-0.5 * (dist / hi)^2)     # Gaussian kernel weights
    num  <- rowSums(W * x[idx])
    den  <- rowSums(W)
    num / den
  }

  ## ---------- bias-variance curves over k ----------
  bw_seq <- 12:max_bandwidth

  LM_A <- sapply(bw_seq, function(k) lw_mean(xA_t, k))
  LM_B <- sapply(bw_seq, function(k) lw_mean(xB_t, k))

  # Bias: mean squared deviation between local mean and original value
  Bias_A <- colMeans((LM_A - xA_t)^2)
  Bias_B <- colMeans((LM_B - xB_t)^2)

  # Variance: variance across polygons of the local mean surface
  Var_A  <- apply(LM_A, 2, var)
  Var_B  <- apply(LM_B, 2, var)

  # standardize to [0,1] for trade-off comparison
  Bias_A_s <- scales::rescale(Bias_A)
  Bias_B_s <- scales::rescale(Bias_B)
  Var_A_s  <- scales::rescale(Var_A)
  Var_B_s  <- scales::rescale(Var_B)

  trade_A <- data.frame(Bandwidth = bw_seq, Bias = Bias_A_s, Variance = Var_A_s)
  trade_B <- data.frame(Bandwidth = bw_seq, Bias = Bias_B_s, Variance = Var_B_s)

  # bandwidth where |Bias - Var| is minimized for each map
  kA <- trade_A$Bandwidth[ which.min(abs(trade_A$Bias - trade_A$Variance)) ]
  kB <- trade_B$Bandwidth[ which.min(abs(trade_B$Bias - trade_B$Variance)) ]

  k_sqrt     <- ceiling(sqrt(n))
  k_rect_min <- min(kA, kB)
  k_rect_max <- max(kA, kB)

  # choose bandwidth
  k_star <-
    if (option == "midpoint") ceiling(mean(c(kA, kB))) else
      if (option == "lower")  k_rect_min else
        if (option == "upper")  k_rect_max else
          stop("option must be one of 'midpoint', 'lower', 'upper'")

  ## ---------- plot (bias-variance trade-off) ----------
  p <- ggplot2::ggplot() +
    # Variable A: Green (Solid for Bias, Dashed for Variance)
    ggplot2::geom_line(data = trade_A,
                       ggplot2::aes(Bandwidth, Bias, color = "Variable A"),
                       linewidth = 0.8) +
    ggplot2::geom_line(data = trade_A,
                       ggplot2::aes(Bandwidth, Variance, color = "Variable A"),
                       linetype = "dashed", linewidth = 0.8) +

    # Variable B: Blue (Solid for Bias, Dashed for Variance)
    ggplot2::geom_line(data = trade_B,
                       ggplot2::aes(Bandwidth, Bias, color = "Variable B"),
                       linewidth = 0.8) +
    ggplot2::geom_line(data = trade_B,
                       ggplot2::aes(Bandwidth, Variance, color = "Variable B"),
                       linetype = "dashed", linewidth = 0.8) +

    # Suggested Range and Reference Lines
    ggplot2::geom_rect(
      ggplot2::aes(xmin = k_rect_min, xmax = k_rect_max, ymin = 0, ymax = 1),
      fill = "grey70", alpha = 0.3
    ) +
    ggplot2::geom_vline(xintercept = k_sqrt, linewidth = 0.5, color = "grey40") +
    ggplot2::geom_vline(xintercept = c(kA, kB), color = "red",
                        linewidth = 0.3, alpha = 0.6) +
    ggplot2::geom_vline(xintercept = k_star, color = "purple",
                        linewidth = 0.8) +

    # Labels and Scales
    ggplot2::geom_text(
      ggplot2::aes(x = k_sqrt, y = 1.02),
      label = paste0("sqrt(n)=", k_sqrt),
      size = 3.5, color = "grey40"
    ) +
    ggplot2::scale_color_manual(
      name = "Variables",
      values = c("Variable A" = "forestgreen", "Variable B" = "royalblue")
    ) +
    ggplot2::labs(
      x = "Bandwidth (k-nearest neighbours)",
      y = "Standardized bias/variance",
      subtitle = paste0("Range [", k_rect_min, ", ", k_rect_max, "]; Chosen = ", k_star),
      caption = "* Solid line: Bias / Dashed line: Variance"
    ) +
    ggplot2::theme_minimal()

  out <- list(
    plot = p,
    bandwidth = k_star,
    tradeoff = list(
      map1   = trade_A,
      map2   = trade_B,
      kA     = kA,
      kB     = kB,
      k_sqrt = k_sqrt
    )
  )

  return(out)
}
utils::globalVariables(c("Bandwidth", "Bias", "Variance"))

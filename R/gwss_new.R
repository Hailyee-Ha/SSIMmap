#' @importFrom sf st_crs st_centroid st_coordinates
#' @importFrom stats cov.wt

gwss_new<-function (data, vars,bw){
  p4s <- st_crs(data)$proj4string
  options(sf_warn_centroid_attrs = FALSE)
  dp.locat <- st_centroid(data)
  dp.locat <- data.frame(st_coordinates(dp.locat))
  dp.locat<-as.matrix(dp.locat)
  sp.locat <- st_centroid(data)
  sp.locat <- data.frame(st_coordinates(sp.locat))
  sp.locat<-as.matrix(sp.locat)

  shape<-data
  data <- as.data.frame(data)
  dp.n <- nrow(data)
  sp.n <- nrow(sp.locat)

  if (missing(vars))
    stop("Variables input error")
  if (missing(bw) || bw <= 0)
    stop("Bandwidth is not specified incorrectly")
  len.var <- length(vars)
  col.nm <- colnames(data)
  var.idx <- match(vars, col.nm)[!is.na(match(vars, col.nm))]
  if (length(var.idx) == 0)
    stop("Variables input doesn't match with data")
  x <- data[, var.idx]
  x <- as.matrix(x)
  var.nms <- names(data)[var.idx]
  var.n <- ncol(x)
  if (len.var > var.n)
    warning("Invalid variables have been specified, please check them again!")
  local.mean <- matrix(numeric(var.n * sp.n), ncol = var.n)
  standard.deviation <- matrix(numeric(var.n * sp.n), ncol = var.n)
  LVar <- matrix(numeric(var.n * sp.n), ncol = var.n)
  cov.nms <- c()
  cov.mat <- c()
  if (var.n > 1) {
    cov.mat <- matrix(numeric((var.n - 1) * var.n * sp.n/2),
                      nrow = sp.n)
  }

  for (i in 1:sp.n) {
    dist.vi <- gw.dist(dp.locat, sp.locat, focus = i)
    W.i <- as.matrix(gw.weight(dist.vi, bw), nrow = 1)
    sum.w <- sum(W.i)
    Wi <- W.i/sum.w
    Wi <- as.numeric(Wi)
    local.mean[i, ] <- Wi %*% x
    for (j in 1:var.n) {
      LVar[i, j] <- Wi %*% ((x[, j] - local.mean[i, j])^2)
      standard.deviation[i, j] <- sqrt(LVar[i, j])
    }
      tag <- 0
      tag <- tag + 1
      cov.mat[i, tag] <- cov.wt(cbind(x[, 1], x[, 2]), wt =Wi)$cov[1, 2]
     }


  colnames(local.mean) <- paste(var.nms, "LM", sep = "_")
  colnames(standard.deviation) <- paste(var.nms, "LSD", sep = "_")
  colnames(LVar) <- paste(var.nms, "LVar", sep = "_")
  if (var.n > 1) {
    for (i in 1:(var.n - 1)) {
      for (j in (i + 1):var.n) {
        cov.v1v2 <- paste("Cov", paste(var.nms[i], var.nms[j],
                                       sep = "."), sep = "_")
        cov.nms <- c(cov.nms, cov.v1v2)
      }
    }
    colnames(cov.mat) <- cov.nms
  }
  res.df <- data.frame(local.mean, standard.deviation,
                           LVar, cov.mat)

  rownames(res.df) <- rownames(sp.locat)

 gwresult<- cbind(shape,res.df)
 return(gwresult)
}

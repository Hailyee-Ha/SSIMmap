#' Calculate SSIM index for raster images
#'
#' This function calculates the SSIM, a measure of similarity between two raster images
#'
#' @param img1 Raster object representing the first image.
#' @param img2 Raster object representing the second image.
#' @param global If global is True, returning the global average of SSIM, SIM, SIV, and SIP. If the option is FALSE, a a \code{raster} raster brick containing the SSIM, SIM, SIV, and SIP for each cell is returned
#' Default is TRUE
#' @param w Integer specifying the window size for the local neighborhood. Default is 3.
#' @param k1 the constant used in the SSIM calculation. Default is NULL, in which case it is computed from the maximum value of variables.
#' @param k2 the constant used in the SSIM calculation. Default is NULL, in which case it is computed from the maximum value of variables.
#' @return If global is True, returning the global average of SSIM, SIM, SIV, and SIP. If the option is FALSE, a \code{raster} raster brick containing the SSIM, SIM, SIV, and SIP for each cell is returned
#' @details This function computes the SSIM index for two raster images.
#'
#' @importFrom terra  global focal rast crop
#' @examples
#'
#' image_path1 <- system.file("ex","groups2nm.tif", package = "SSIMmap", mustWork = TRUE)
#' image_path2 <- system.file("ex","single2nm.tif", package = "SSIMmap", mustWork = TRUE)
#' img1<-terra::rast(image_path1)
#' img2<-terra::rast(image_path2)
#' ssim_raster(img1,img2)
#' result_raster<-ssim_raster(img1,img2,global=FALSE)
#' @export ssim_raster


ssim_raster<- function(img1, img2, global=TRUE, w=3,k1=NULL,k2=NULL) {
  img1.extent <- terra::ext(img1)
  img1.na <- is.na(img1)

  #set constants
  l <- max(terra::global(img1, fun= max,na.rm=TRUE), terra::global(img2, fun= max, na.rm=TRUE))
  globalMin <- abs(min(terra::global(img1, fun= min,na.rm=TRUE), terra::global(img2, fun= min,na.rm=TRUE)))
  l <- l - globalMin

  if(is.null(k1)){
    k1<-(0.01* l )/255
  }
  if(is.null(k2)){
    k2<-(0.03* l )/255
  }

  if(k1>1){
    stop("k1 needs to be smaller than 1")
  }

  if(k2>1){
    stop("k2 needs to be smaller than 1")
  }
  k <- c(k1, k2)
  C1 <-(k[1]*l)^2
  C2 <-(k[2]*l)^2
  C3 <-C2/2

  #Create Null filter
  filterx <- matrix(1,ncol=w*2+1,nrow=w*2+1)/(w*2+1)^2


  #get mu
  mu1 <- terra::focal(img1, filterx)
  mu2 <- terra::focal(img2, filterx)

  sig1 <- abs(terra::focal(img1*img1,filterx) - mu1*mu1)^0.5
  sig2 <- abs(terra::focal(img2*img2,filterx) - mu2*mu2)^0.5
  #sig12 relates to correlation
  sig12 <- terra::focal(img1*img2, filterx) - mu1*mu2

  #compute components
  L <- ((2*mu1*mu2)+C1) / (mu1^2 + mu2^2 + C1)
  C <- ((2*sig1*sig2)+C2) / (sig1^2 + sig2^2 + C2)
  S <- (sig12 + C3) / (sig1 * sig2 + C3)
  #compute SSIM
  SSIM2 <- L * C * S

  #Compute RasterBrick
  ssim.brick <- c(SSIM2,L,C,S)
  ssim.brick <- terra::crop(ssim.brick, img1.extent)
  ssim.brick[img1.na] <- NA

  names(ssim.brick) <- c("SSIM", "SIM", "SIV", "SIP")
  mean_SSIM <- terra::global(SSIM2, fun= 'mean', na.rm = TRUE)
  mean_SIM <- terra::global(L, fun= 'mean', na.rm = TRUE)
  mean_SIV <- terra::global(C, fun= 'mean', na.rm = TRUE)
  mean_SIP <- terra::global(S, fun= 'mean', na.rm = TRUE)

  if(global){
    result <- paste("SSIM:", round(mean_SSIM, 5),
                    "SIM:", round(mean_SIM, 5),
                    "SIV:", round(mean_SIV, 5),
                    "SIP:", round(mean_SIP, 5),
                    sep = " ")


    cat(result, "\n")}

  else{
    return(ssim.brick)
  }
}


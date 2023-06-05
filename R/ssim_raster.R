ssim_raster<- function(img1, img2, w=3,k1=NULL,k2=NULL) {

  #Check to see if extents are equal
  img1.extent <- extent(img1)
  img2.extent <- extent(img2)
  img1.na <- Which(is.na(img1),cells=TRUE)
  if (img1.extent != img2.extent){stop('Warning: SSIM calculation aborted. The raster extents do not match.')}

  #set constants
  l <- max(cellStats(img1, max), cellStats(img2, max))
  globalMin <- abs(min(cellStats(img1, min), cellStats(img2, min)))
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
  mu1 <- focal(img1, filterx)
  mu2 <- focal(img2, filterx)

  sig1 <- abs(focal(img1*img1,filterx) - mu1*mu1)^0.5
  sig2 <- abs(focal(img2*img2,filterx) - mu2*mu2)^0.5
  #sig12 relates to correlation
  sig12 <- focal(img1*img2, filterx) - mu1*mu2

  #compute components
  L <- ((2*mu1*mu2)+C1) / (mu1^2 + mu2^2 + C1)
  C <- ((2*sig1*sig2)+C2) / (sig1^2 + sig2^2 + C2)
  S <- (sig12 + C3) / (sig1 * sig2 + C3)
  #compute SSIM
  SSIM2 <- L * C * S

  #Compute RasterBrick
  ssim.brick <- brick(SSIM2, L, C, S)
  ssim.brick <- crop(ssim.brick,img1.extent)
  ssim.brick[img1.na] <- NA

  ssim.brick@data@names <- c('SSIM', 'SIM', 'SIV', 'SIP')
  return(ssim.brick)
}


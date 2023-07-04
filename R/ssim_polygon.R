#' Calculate the structural similarity index measure for polygon maps.
#'
#' This function computes the SSIM, a measure of similarity between two polygon maps.
#' The ssim_polygon function computes the SSIM for each polygon and can be returned as a global average or for each polygon as a \code{sf} object
#'
#' @param shape a a \code{sf} polygon containing the polygon data with attributes that can create polygon-based maps
#' @param map1 the name of the first map to compare as a column in the shape
#' @param map2 the name of the second map to compare as a column in the shape
#' @param global If global is True, returning the global average of SSIM, SIM, SIV, and SIP. If the option is FALSE, a a \code{sf} SpatialPolygonsDataFrame containing the SSIM, SIM, SIV, and SIP for each polygon is returned
#' Default is TRUE
#' @param k1 the constant used in the SSIM calculation. Default is NULL, in which case it is computed from the maximum value of variables.
#' @param k2 the constant used in the SSIM calculation. Default is NULL, in which case it is computed from the maximum value of variables.
#' @param bandwidth bandwidth for the Gaussian kernel weighting used in the SSIM calculation. Default is the square root of N
#' @param standardize If TRUE, standardize the variables before computing the SSIM. Default is TRUE.
#'
#'
#' @details This function computes the SSIM index for two polygon maps.
#' @return If global is TRUE, a string containing the global average SSIM, SIM, SIV, and SIP.
#' If global is FALSE, a \code{sf} polygon containing the SSIM, SIM, SIV, and SIP for each polygon.
#'
#' @importFrom sf as_Spatial
#' @importFrom dplyr select contains
#' @importFrom GWmodel gwss
#' @importFrom knitr kable

#' @examples
#' # Load example sf polygon Toronto Areas with attributes for maps:
#' # Pampalon Index,CIMD Index,
#' #and percentage of household commuting within the same Census Sub Division of residence
#' shape<-SSIMmap::polygon
#'
#' # Mapping two attributes
#' plot(shape$CIMD_SDD)
#' plot(shape$PP_SDD)
#' #Finding global ssim
#'
#' ssim_polygon(shape,"CIMD_SDD","PP_SDD")
#' #Finding local ssim
#' df<-ssim_polygon(shape,"CIMD_SDD","PP_SDD",global=FALSE)

#' @export ssim_polygon


# Functions ---------------------------------------------------------------

ssim_polygon<-function(shape,map1,map2,standardize=TRUE,bandwidth=NULL,k1=NULL,k2=NULL,global=TRUE){
  if(map1==map2){
    stop("variables are identical")
  }

  num_rows <- nrow(shape)
  sqrt_num_rows <- round(sqrt(num_rows),0)


  if(standardize){
    shape_df<-as.data.frame(shape)
    shape_df$z_score_map1<-(shape_df[,map1]-mean(shape_df[,map1]))/sd(shape_df[,map1])
    shape_df$z_score_map2<-(shape_df[,map2]-mean(shape_df[,map2]))/sd(shape_df[,map2])
    shape_df$z_score_map1<-shape_df$z_score_map1-min(shape_df$z_score_map1)
    shape_df$z_score_map2<-shape_df$z_score_map2-min(shape_df$z_score_map2)
    z_scores<-as.data.frame(cbind(shape_df$z_score_map1,shape_df$z_score_map2))
    colnames(z_scores)<- paste0("zscore",colnames(z_scores))
    names(z_scores)<- sub("V1",map1,names(z_scores))
    names(z_scores)<- sub("V2",map2,names(z_scores))
    shape_merged<-cbind(shape,z_scores)
    shape_merged<-sf::as_Spatial(shape_merged)
    map1<-as.character(colnames(z_scores[1]))
    map2<-as.character(colnames(z_scores[2]))
    if(is.null(bandwidth)){
      result<-GWmodel::gwss(shape_merged,vars = c(map1,map2), kernel = "gaussian",adaptive = TRUE,bw=sqrt_num_rows)

    }
    else{
      result<-GWmodel::gwss(shape_merged,vars = c(map1,map2), kernel = "gaussian",adaptive = TRUE,bw=bandwidth)
    }
    gwss_result<-as.data.frame(result$SDF)
    mean<-dplyr::select(gwss_result,contains("LM"))
    sd<-dplyr::select(gwss_result,contains("LSD"))
    cov<-dplyr::select(gwss_result,contains("Cov"))
    merged_df<-cbind(mean,sd,cov,z_scores)
    merged_df<-cbind(merged_df,shape)



  }

  else{
    shape_sp<-sf::as_Spatial(shape)

    if(is.null(bandwidth)){
      result<-GWmodel::gwss(shape_merged,vars = c(map1,map2), kernel = "gaussian",adaptive = TRUE,bw=sqrt_num_rows)

    }else{
    result<-GWmodel::gwss(shape_sp,vars = c(map1,map2), kernel = "gaussian",adaptive = TRUE,bw=bandwidth)
    }
    gwss_result<-as.data.frame(result$SDF)
    mean<-dplyr::select(gwss_result,contains("LM"))
    sd<-dplyr::select(gwss_result,contains("LSD"))
    cov<-dplyr::select(gwss_result,contains("Cov"))
    merged_df<-cbind(mean,sd,cov)
    colnames(merged_df)<-gsub("zscore", "", colnames(merged_df))
    merged_df<-cbind(merged_df,shape)




  }


  shape_df<-as.data.frame(merged_df)
  map1<-as.data.frame(dplyr::select(shape_df,contains(map1)))
  map2<-as.data.frame(dplyr::select(shape_df,contains(map2)))

  globalMin <- min(map1,map2)
  l<-max(map1,map2)
  max<-max(map1,map2)

  map1_LM<-dplyr::select(map1,contains('_LM'))
  map2_LM<-dplyr::select(map2,contains('_LM'))
  map1_LSD<-dplyr::select(map1,contains('_LSD'))
  map2_LSD<-dplyr::select(map2,contains('_LSD'))
  Cov_map1_map2<-dplyr::select(shape_df,contains('Cov_'))
  l <- l-globalMin
  options(scipen=999)


  if(is.null(k1)){
    k1<-(0.01*max)/255
  }
  if(is.null(k2)){
    k2<-(0.03*max)/255
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

  mu1<-map1_LM
  mu2<-map2_LM
  sig1<-map1_LSD
  sig2<-map2_LSD
  sig12<-Cov_map1_map2
  SIM<- ((2*mu1*mu2)+C1) / (mu1^2 + mu2^2 + C1)
  if(any(SIM>1) || any(SIM<0)){
    stop("violating the boundness rule,you need to rescale the attributes for maps")
  }

  SIV <- ((2*sig1*sig2)+C2) / (sig1^2 + sig2^2 + C2)
  if(any(SIV>1) || any(SIV<0)){
    stop("violating the boundness rule,you need to rescale the attributes for maps")
  }
  SIP <- (sig12 + C3) / (sig1 * sig2 + C3)
  if(any(SIP>1)|| any(SIP< (-1))){
    stop("violating the boundness rule,you need to rescale the attributes for maps")
  }
  SSIM <- SIM*SIV*SIP
  df<-as.data.frame(cbind(SSIM, SIM, SIV, SIP))
  colnames(df)<-c("SSIM","SIM","SIV","SIP")

  if(global){
    SSIM_mean <- mean(df$SSIM)
    SIM_mean <- mean(df$SIM)
    SIV_mean <- mean(df$SIV)
    SIP_mean <- mean(df$SIP)

    SSIM_min <- min(df$SSIM)
    SIM_min <- min(df$SIM)
    SIV_min <- min(df$SIV)
    SIP_min <- min(df$SIP)

    SSIM_max <- max(df$SSIM)
    SIM_max <- max(df$SIM)
    SIV_max <- max(df$SIV)
    SIP_max <- max(df$SIP)

    SSIM_sd <- sd(df$SSIM)
    SIM_sd <- sd(df$SIM)
    SIV_sd <- sd(df$SIV)
    SIP_sd <- sd(df$SIP)

    result <- data.frame(
      Statistic = c("Mean", "Min", "Max", "SD"),
      SSIM = c(SSIM_mean, SSIM_min, SSIM_max, SSIM_sd),
      SIM = c(SIM_mean, SIM_min, SIM_max, SIM_sd),
      SIV = c(SIV_mean, SIV_min, SIV_max, SIV_sd),
      SIP = c(SIP_mean, SIP_min, SIP_max, SIP_sd)
    )

    print(knitr::kable(result, row.names = FALSE))
  }
  else{
    SSIM_results<-cbind(shape,df)
    return(SSIM_results)
  }
}

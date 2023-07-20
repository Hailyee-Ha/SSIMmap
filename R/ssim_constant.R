#' Calculate constants for SSIM Index
#'
#' This function rescales the constants (k1 and k2) for the SSIM index based on the global maximum value of the maps
#' @param shape a a \code{sf} polygon  containing the polygon data with attributes that can create polygon-based maps
#' @param map1 the name of the first map to compare
#' @param map2 the name of the second map to compare
#' @param standardize If TRUE, standardize the variables before computing the SSIM. Default is TRUE.
#' @importFrom stats sd var mean
#' @return The rescaled constants (k1 and k2)
#' @details This function calculates the rescaled constants (k1 and k2) for the SSIM index.
#' k1 and k2 in the original SSIM method,which are for the 8-bit grey scale images, are 0.01 and 0.03 respectively.
#' The SSIM for maps can use the rescaled k1 and k2 based on the global maximum value of two maps.
#'
#' @examples
#' # Load example sf class Toronto Area with attributes for maps:
#' # Toronto Areas with attributes for maps:Pampalon Index,CIMD Index,
#' # and percentage of household commuting within the same Census Sub Division of residence)
#' shape<-SSIMmap::Toronto
#' ssim_constant(shape,"PP_SDD","CIMD_SDD")
#'
#' @export ssim_constant

# Functions ---------------------------------------------------------------

ssim_constant<-function(shape,map1,map2,standardize=TRUE){
  if(map1==map2){
    stop("maps are identical")
  }

  if(standardize){
    shape_df<-as.data.frame(shape)
    shape_df$z_score_map1<-(shape_df[,map1]-mean(shape_df[,map1]))/sd(shape_df[,map1])
    shape_df$z_score_map2<-(shape_df[,map2]-mean(shape_df[,map2]))/sd(shape_df[,map2])
    min<-min(min(shape_df$z_score_map1),min(shape_df$z_score_map2))
    shape_df$z_score_map1<-shape_df$z_score_map1-min
    shape_df$z_score_map2<-shape_df$z_score_map2-min
    max<-max(shape_df$z_score_map1,shape_df$z_score_map2)
  }

  else{
    shape_df<-as.data.frame(shape)
    max<-max(shape_df[,map1],shape_df[,map2])
  }


  k1<-round((0.01*max)/255,5)
  k2<-round((0.03*max)/255,5)
  result<- paste("Rescaled K1:", k1,"Rescaled K2:", k2 ,sep = " ")
  cat(result)
}


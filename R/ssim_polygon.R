ssim_polygon<-function(shape,var1,var2,global=TRUE,k1=NULL,k2=NULL,bandwidth=NULL,standardize=TRUE){
  if(var1==var2){
    stop("variables are identical")
  }
  if(bandwidth < 12){
    stop("The Bandwdith is too small to calculate the SSIM for the irregular lattice maps")
  }

  if(standardize){
    shape_df<-as.data.frame(shape)
    shape_df$z_score_var1<-(shape_df[,var1]-mean(shape_df[,var1]))/sd(shape_df[,var1])
    shape_df$z_score_var2<-(shape_df[,var2]-mean(shape_df[,var2]))/sd(shape_df[,var2])
    min<-min(min(shape_df$z_score_var1),min(shape_df$z_score_var2))
    shape_df$z_score_var1<-shape_df$z_score_var1-min
    shape_df$z_score_var2<-shape_df$z_score_var2-min
    z_scores<-as.data.frame(cbind(shape_df$z_score_var1,shape_df$z_score_var2))
    colnames(z_scores)<- paste0("zscore",colnames(z_scores))
    names(z_scores)<- sub("V1",var1,names(z_scores))
    names(z_scores)<- sub("V2",var2,names(z_scores))
    shape_merged<-cbind(shape,z_scores)
    shape_merged<-as(shape_merged,"Spatial")
    var1<-as.character(colnames(z_scores[1]))
    var2<-as.character(colnames(z_scores[2]))
    result<-GWmodel::gwss(shape_merged,vars = c(var1,var2), kernel = "gaussian",adaptive = TRUE,bw=bandwidth)
    gwss_result<-as.data.frame(result$SDF)
    mean<-dplyr::select(gwss_result,contains("LM"))
    sd<-dplyr::select(gwss_result,contains("LSD"))
    cov<-dplyr::select(gwss_result,contains("Cov"))
    merged_df<-cbind(mean,sd,cov,z_scores)
    merged_df<-cbind(merged_df,shape)



  }

  else{
    shape_sp<-as(shape,"Spatial")
    result<-GWmodel::gwss(shape_sp,vars = c(var1,var2), kernel = "gaussian",adaptive = TRUE,bw=bandwidth)
    gwss_result<-as.data.frame(result$SDF)
    mean<-dplyr::select(gwss_result,contains("LM"))
    sd<-dplyr::select(gwss_result,contains("LSD"))
    cov<-dplyr::select(gwss_result,contains("Cov"))
    merged_df<-cbind(mean,sd,cov)
    colnames(merged_df)<-gsub("zscore", "", colnames(merged_df))
    merged_df<-cbind(merged_df,shape)




  }


  shape_df<-as.data.frame(merged_df)
  var1<-as.data.frame(dplyr::select(shape_df,contains(var1)))
  var2<-as.data.frame(dplyr::select(shape_df,contains(var2)))

  globalMin <- min(var1,var2)
  l<-max(var1,var2)
  max<-max(var1,var2)

  var1_LM<-dplyr::select(var1,contains('_LM'))
  var2_LM<-dplyr::select(var2,contains('_LM'))
  var1_LSD<-dplyr::select(var1,contains('_LSD'))
  var2_LSD<-dplyr::select(var2,contains('_LSD'))
  Cov_var1_var2<-dplyr::select(shape_df,contains('Cov_'))
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

  mu1<-var1_LM
  mu2<-var2_LM
  sig1<-var1_LSD
  sig2<-var2_LSD
  sig12<-Cov_var1_var2
  SIM<- ((2*mu1*mu2)+C1) / (mu1^2 + mu2^2 + C1)
  if(any(SIM>1) || any(SIM<0)){
    stop("violating the boundness rule,you need to rescale the variables")
  }

  SIV <- ((2*sig1*sig2)+C2) / (sig1^2 + sig2^2 + C2)
  if(any(SIV>1) || any(SIV<0)){
    stop("violating the boundness rule,you need to rescale the variables")
  }
  SIP <- (sig12 + C3) / (sig1 * sig2 + C3)
  if(any(SIP>1)|| any(SIP< (-1))){
    stop("violating the boundness rule,you need to rescale the variables")
  }
  SSIM <- SIM*SIV*SIP
  df<-as.data.frame(cbind(SSIM, SIM, SIV, SIP))
  colnames(df)<-c("SSIM","SIM","SIV","SIP")

  if(global){
    SSIM<-mean(df$SSIM)
    SIM<-mean(df$SIM)
    SIV<-mean(df$SIV)
    SIP<-mean(df$SIP)
    result <- paste("SSIM:", round(SSIM, 5),
                    "SIM:", round(SIM, 5),
                    "SIV:", round(SIV, 5),
                    "SIP:", round(SIP, 5),
                    sep = " ")


    cat(result, "\n")}

  else{
    SSIM_results<-cbind(shape,df)
    return(SSIM_results)
  }
}

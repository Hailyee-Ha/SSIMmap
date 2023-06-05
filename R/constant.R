constant<-function(shape,var1,var2,standardize=TRUE){
  if(var1==var2){
    stop("variables are identical")
  }

  if(standardize){
    shape_df<-as.data.frame(shape)
    shape_df$z_score_var1<-(shape_df[,var1]-mean(shape_df[,var1]))/sd(shape_df[,var1])
    shape_df$z_score_var2<-(shape_df[,var2]-mean(shape_df[,var2]))/sd(shape_df[,var2])
    min<-min(min(shape_df$z_score_var1),min(shape_df$z_score_var2))
    shape_df$z_score_var1<-shape_df$z_score_var1-min
    shape_df$z_score_var2<-shape_df$z_score_var2-min
    max<-max(shape_df$z_score_var1,shape_df$z_score_var2)
  }

  else{
    shape_df<-as.data.frame(shape)
    max<-max(shape_df[,var1],shape_df[,var2])
  }


  k1<-round((0.01*max)/255,5)
  k2<-round((0.03*max)/255,5)
  result<- paste("Rescaled K1:", k1,"Rescaled K2:", k2 ,sep = " ")
  cat(result)
}


bandwidth<-function(shape, var1,var2,max_bandwidth=max_bandwidth,standarize=TRUE,option="median"){

  if(var1==var2){
    stop("variables are identical")
  }
  if(max_bandwidth < 12){
    stop("The maximum size of the bandwdith is too small to compute the SSIM index")
  }
  if(nrow(shape)<max_bandwidth){
    stop("The maximum size of the bandwdith is larger than the number of lattice")
  }

  if(standarize){
    shape_df<-as.data.frame(shape)
    shape_df$z_score_var1<-(shape_df[,var1]-mean(shape_df[,var1]))/sd(shape_df[,var1])
    shape_df$z_score_var2<-(shape_df[,var2]-mean(shape_df[,var2]))/sd(shape_df[,var2])
    min<-min(min(shape_df$z_score_var1),min(shape_df$z_score_var2))
    shape_df$z_score_var1<-shape_df$z_score_var1-min
    shape_df$z_score_var2<-shape_df$z_score_var2-min
    mean_shape_var1<-mean(shape_df$z_score_var1)
    mean_shape_var2<-mean(shape_df$z_score_var2)
    var_shape_var1<-sd(shape_df$z_score_var1)
    var_shape_var2<-sd(shape_df$z_score_var2)
    cov_shape<-cov(shape_df$z_score_var1,shape_df$z_score_var2)
    z_scores<-as.data.frame(cbind(shape_df$z_score_var1,shape_df$z_score_var2))
    colnames(z_scores)<- paste0("zscore",colnames(z_scores))
    names(z_scores)<- sub("V1",var1,names(z_scores))
    names(z_scores)<- sub("V2",var2,names(z_scores))
    var1<-as.character(colnames(z_scores[1]))
    var2<-as.character(colnames(z_scores[2]))
    globalMin <- min(min(z_scores[,1]),min(z_scores[,2]))
    l<-max(max(z_scores[,1]),max(z_scores[,2]))
    max<-max(max(z_scores[,1]),max(z_scores[,2]))
    shape_merged<-cbind(shape,z_scores)


    bw=seq(from=12,to=max_bandwidth,by=1)
    shape_merged<-as(shape_merged,"Spatial")
    gwsses<-list()

    for(num in bw){
      result<-GWmodel::gwss(shape_merged,vars = c(var1,var2), kernel = "gaussian",adaptive = TRUE,bw=num)
      gwsses[[num]]<-result$SDF
    }
    df_t<-as.data.frame(matrix())

    for(i in bw){
      temp_df<-as.data.frame(gwsses[[i]])
      colnames(temp_df)<-c(paste(colnames(temp_df),i,sep="."))
      df_t<-cbind(df_t,temp_df)
    }
    df_t<-df_t[-1]

    var1_o<-shape_df$z_score_var1
    var2_o<-shape_df$z_score_var2
  }
  else{
    shape_df<-as.data.frame(shape)
    bw=seq(from=12,to=max_bandwidth,by=1)
    gwsses<-list()
    for(num in bw){
      result<-GWmodel::gwss(shape,vars = c(var1,var2), kernel = "gaussian",adaptive = TRUE,bw=num)
      gwsses[[num]]<-result$SDF
    }
    df_t<-as.data.frame(matrix())

    for(i in bw){
      temp_df<-as.data.frame(gwsses[[i]])
      colnames(temp_df)<-c(paste(colnames(temp_df),i,sep="."))
      df_t<-cbind(df_t,temp_df)
    }
    df_t<-df_t[-1]
    var1_o<-shape_df$var1
    var2_o<-shape_df$var2
  }
  var1<-as.data.frame(dplyr::select(df_t,contains(var1)))
  var2<-as.data.frame(dplyr::select(df_t,contains(var2)))
  var1_LM<-dplyr::select(var1,contains('_LM'))
  var2_LM<-dplyr::select(var2,contains('_LM'))

  col_names<-paste(seq(from=12,to=max_bandwidth,by=1))
  colnames(var1_LM)<-col_names
  colnames(var2_LM)<-col_names
  df_bias_var1<-apply(var1_LM,2,function(x)(x-var1_o)^2)
  df_bias_var1<-as.data.frame(apply(df_bias_var1,2, mean))
  df_bias_var2<-apply(var2_LM,2,function(x)(x-var2_o)^2)
  df_bias_var2<-as.data.frame(apply(df_bias_var2,2, mean))
  df_variance_var1<-as.data.frame(apply(var1_LM,2,var))
  df_variance_var2<-as.data.frame(apply(var2_LM,2,var))
  bw_order<-as.data.frame(bw)

  df_bias_var1<-as.data.frame(cbind(bw_order,df_bias_var1))
  df_bias_var2<-as.data.frame(cbind(bw_order,df_bias_var2))
  df_variance_var1<-as.data.frame(cbind(bw_order,df_variance_var1))
  df_variance_var2<-as.data.frame(cbind(bw_order,df_variance_var2))
  colnames(df_bias_var1)<-c("Bandwidth","Bias")
  colnames(df_bias_var2)<-c("Bandwidth","Bias")
  colnames(df_variance_var1)<-c("Bandwidth","Variance")
  colnames(df_variance_var2)<-c("Bandwidth","Variance")
  df_bias_var1$R_Bias<-rescale(df_bias_var1$Bias)
  df_bias_var2$R_Bias<-rescale(df_bias_var2$Bias)
  df_variance_var1$R_Variance<-rescale(df_variance_var1$Variance)
  df_variance_var2$R_Variance<-rescale(df_variance_var2$Variance)

  P_Tradeoff_var1<-as.data.frame(df_bias_var1$R_Bias-df_variance_var1$R_Variance)
  Tradeoff_var1<-as.data.frame(cbind(bw_order,P_Tradeoff_var1))
  colnames(Tradeoff_var1)<-c("Bandwidth","Tradeoff")
  index_var1<-which.min(abs(Tradeoff_var1$Tradeoff))
  bw_closest_to_zero_var1<-Tradeoff_var1$Bandwidth[index_var1]

  P_Tradeoff_var2<-as.data.frame(df_bias_var2$R_Bias-df_variance_var2$R_Variance)
  Tradeoff_var2<-as.data.frame(cbind(bw_order,P_Tradeoff_var2))
  colnames(Tradeoff_var2)<-c("Bandwidth","Tradeoff")
  index_var2<-which.min(abs(Tradeoff_var2$Tradeoff))
  bw_closest_to_zero_var2<-Tradeoff_var2$Bandwidth[index_var2]

  num_rows <- nrow(shape)
  sqrt_num_rows <- round(sqrt(num_rows),1)

  plot<-ggplot()+geom_line(data=df_bias_var1,aes(Bandwidth,R_Bias),color="dark blue")+geom_line(data=df_variance_var1,aes(Bandwidth,R_Variance),color="dark green")+geom_line(data=df_bias_var2,aes(Bandwidth,R_Bias),linetype="dashed",color="dark blue")+geom_line(data=df_variance_var2,aes(Bandwidth,R_Variance),linetype="dashed",color="dark green")
  plot<-plot+geom_rect(aes(xmin = min(bw_closest_to_zero_var1,bw_closest_to_zero_var2), xmax = max(bw_closest_to_zero_var1,bw_closest_to_zero_var2), ymin = 0, ymax = 1), fill = "grey", alpha = 0.5)
  plot<-plot+geom_vline(xintercept = sqrt_num_rows, color="black",linewidth=0.3, alpha = 0.5)+geom_vline(xintercept = bw_closest_to_zero_var1,color="red",linewidth=0.3, alpha = 0.5)+geom_vline(xintercept = bw_closest_to_zero_var2,color="red",linewidth=0.3, alpha = 0.5)+labs(x="Bandwidth",y="bias/variance")
  plot+geom_text(aes(x= bw_closest_to_zero_var1,y=0.5), label = as.character(bw_closest_to_zero_var1),vjust = -1)+geom_text(aes(x= bw_closest_to_zero_var2,y=0.5), label = as.character(bw_closest_to_zero_var2),vjust = -1)+geom_text(aes(x= sqrt_num_rows,y=0.3), label = as.character(sqrt_num_rows),vjust = -1)

  if(option=="median"){
    result<- paste("Square root N:", sqrt_num_rows,"Bias-Variance Trade-off:",mean(c(bw_closest_to_zero_var1,bw_closest_to_zero_var2)),sep = " ")
  }

  else if(option=="lower"){
    result<- paste("Square root N:", sqrt_num_rows,"Bias-Variance Trade-off:",min(c(bw_closest_to_zero_var1,bw_closest_to_zero_var2)),sep = " ")
  }

  else if(option=="upper"){
    result<- paste("Square root N:", sqrt_num_rows,"Bias-Variance Trade-off:",max(c(bw_closest_to_zero_var1,bw_closest_to_zero_var2)),sep = " ")
  }

  else{
    result<-"Invalid option"
  }
  cat(result)
  return(plot)

}



gw.weight<- function(vdist,bw){

  if (is.matrix(vdist)){
    nr <- nrow(vdist)
    dn <- bw/nr
    if(dn<=1)
    {
      rnk<-apply(vdist,2,rank,ties.method='first')
      bw<- vdist[rnk==bw]
    }
    else
    {
      bw <- dn*apply(vdist,2,max)
    }
    if(length(bw)>0)
      wgt<- exp(vdist^2 / (-2 * bw^2))
    else
      wgt <- diag(1, dim(vdist)[1], dim(vdist)[2])
  }else{
    nr <- length(vdist)
    dn <- bw/nr
    if(dn<=1)
    {
      rnk<-rank(vdist,ties.method='first')
      cond<- which(rnk <= bw)
      bw<- vdist[rnk==bw]
    }
    else
    {
      bw <- dn*max(vdist)
    }
      wgt<- exp(vdist^2 / (-2 * bw^2))
  wgt
  }
}


#==================================================#
#   ESTIMATED CV ERROR FOR ANNUAL MODEL            #
#==================================================#
fn_cverr_annu<-function(DATA,POST) {
  # INPUT   DATA (data.frame): data for one measure
  #         POST (list): posterior draws for one measure
  # OUTPUT  int: log pseudomarginal likelihood
  years<-names(table(DATA$year))
  dlen<-ncol(POST[[1]]$eta)
  pyy<-rep(NA,ncol(POST[[1]]$eta)*length(years))
  for (y in 1:length(POST)) {
    yDATA<-DATA[DATA$year==years[y],]
    yPOST<-POST[[y]]
    R<-nrow(yPOST$eta)
    g<-length(yPOST$theta)
    for (d in 1:ncol(yPOST$eta)) {
      mu2<-yPOST$eta[,d]-yPOST$theta[[g]][,d]
      w<-1/(dnorm(yDATA$Ybar[d],mean=mu2,sd=sqrt(yDATA$V[d]+yPOST$sigsq[[g]])))
      pyy[dlen*(y-1)+d]<-R/sum(w)
    }
  }
  return(sum(log(pyy)))
}

#==================================================#
#   ESTIMATED CV ERROR FOR ALL OTHER MODELS        #
#==================================================#
fn_cverr<-function(DATA,POST) {
  # INPUT   DATA (data.frame): data for one measure
  #         POST (list): posterior draws for one measure
  # OUTPUT  int: log pseudomarginal likelihood
  R<-nrow(POST$eta)
  pyy<-rep(NA,ncol(POST$eta))
  g<-length(POST$theta)
  for (d in 1:ncol(POST$eta)) {
    mu2<-POST$eta[,d]-POST$theta[[g]][,d]
    w<-1/(dnorm(DATA$Ybar[d],mean=mu2,sd=sqrt(DATA$V[d]+POST$sigsq[[g]])))
    pyy[d]<-R/sum(w)
  }
  return(sum(log(pyy)))
}

#==================================================#
#   POSTERIOR PREDICTIVE CHECKS                    #
#==================================================#
fn_ppc<-function(ETA,DATA) {
  # INPUT   ETA (vector): posterior draws of eta
  #         DATA (data.frame): data for one measure
  # OUTPUT  matrix: predictive draws
  Yrep<-matrix(NA,nrow=nrow(ETA),ncol=ncol(ETA))
  for (r in 1:nrow(ETA)) {
    Yrep[r,]<-rnorm(ncol(ETA),mean=ETA[r,],sd=sqrt(DATA$V))
  }
  return(Yrep)
}

#==================================================#
#   PERFORMING MODEL COMPARISON                    #
#==================================================#
fn_model_comparison<-function() {
  # OUTPUT   CVerrs (data.frame): log pseudomarginal likelihood for the models
  #          PV (data.frame): posterior predictive p-values
  measures<-fn_data_prep()
  
  annu<-list()
  base<-list()
  sby<-list()
  rich<-list()
  full<-list()
  
  for (m in names(measures)) {
    annu[[m]]<-readRDS(paste('posterior/annu_post_',m,'.rds',sep=''))
    base[[m]]<-readRDS(paste('posterior/base_post_',m,'.rds',sep=''))
    sby[[m]]<-readRDS(paste('posterior/sby_post_',m,'.rds',sep=''))
    rich[[m]]<-readRDS(paste('posterior/rich_post_',m,'.rds',sep=''))
    full[[m]]<-readRDS(paste('posterior/full_post_',m,'.rds',sep=''))
  }
  
  CVerrs<-matrix(NA,nrow=length(measures),ncol=5)
  colnames(CVerrs)<-c('annu','base','sby','rich','full')
  rownames(CVerrs)<-names(measures)
  
  for (m in 1:length(measures)) {
    CVerrs[m,1]<-fn_cverr_annu(DATA=measures[[names(measures)[m]]],POST=annu[[names(annu)[m]]])
    CVerrs[m,2]<-fn_cverr(DATA=measures[[names(measures)[m]]],POST=base[[names(base)[m]]])
    CVerrs[m,3]<-fn_cverr(DATA=measures[[names(measures)[m]]],POST=sby[[names(sby)[m]]])
    CVerrs[m,4]<-fn_cverr(DATA=measures[[names(measures)[m]]],POST=rich[[names(rich)[m]]])
    CVerrs[m,5]<-fn_cverr(DATA=measures[[names(measures)[m]]],POST=full[[names(full)[m]]])
  }
  
  # p-value calculation
  R<-1000
  n<-470
  
  Yannu<-array(NA,dim=c(R,n,length(annu)))
  Ybase<-array(NA,dim=c(R,n,length(base)))
  Ysby<-array(NA,dim=c(R,n,length(sby)))
  Yrich<-array(NA,dim=c(R,n,length(rich)))
  Yfull<-array(NA,dim=c(R,n,length(full)))
  
  for (m in 1:length(base)) {
    annueta<-cbind(annu[[names(annu)[m]]][[1]]$eta,
                   annu[[names(annu)[m]]][[2]]$eta,
                   annu[[names(annu)[m]]][[3]]$eta,
                   annu[[names(annu)[m]]][[4]]$eta,
                   annu[[names(annu)[m]]][[5]]$eta)
    Yannu[,,m]<-fn_ppc(ETA=annueta,DATA=measures[[names(measures)[m]]])
    Ybase[,,m]<-fn_ppc(ETA=base[[names(base)[m]]]$eta,DATA=measures[[names(measures)[m]]])
    Ysby[,,m]<-fn_ppc(ETA=sby[[names(sby)[m]]]$eta,DATA=measures[[names(measures)[m]]])
    Yrich[,,m]<-fn_ppc(ETA=rich[[names(rich)[m]]]$eta,DATA=measures[[names(measures)[m]]])
    Yfull[,,m]<-fn_ppc(ETA=full[[names(full)[m]]]$eta,DATA=measures[[names(measures)[m]]])
  }
  
  PV<-matrix(NA,nrow=length(base),ncol=5)
  
  for (m in 1:length(base)) {
    pv<-round(mean(sd(measures[[names(measures)[m]]]$Ybar)<apply(Yannu[,,m],1,sd)),3)
    PV[m,1]<-min(pv,1-pv)
    pv<-round(mean(sd(measures[[names(measures)[m]]]$Ybar)<apply(Ybase[,,m],1,sd)),3)
    PV[m,2]<-min(pv,1-pv)
    pv<-round(mean(sd(measures[[names(measures)[m]]]$Ybar)<apply(Ysby[,,m],1,sd)),3)
    PV[m,3]<-min(pv,1-pv)
    pv<-round(mean(sd(measures[[names(measures)[m]]]$Ybar)<apply(Yrich[,,m],1,sd)),3)
    PV[m,4]<-min(pv,1-pv)
    pv<-round(mean(sd(measures[[names(measures)[m]]]$Ybar)<apply(Yfull[,,m],1,sd)),3)
    PV[m,5]<-min(pv,1-pv)
  }
  
  colnames(PV)<-c('annu','base','sby','rich','full')
  rownames(PV)<-names(measures)
  
  return(list(CVerrs=CVerrs,PV=PV))
}
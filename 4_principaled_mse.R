#==================================================#
#   DRAWS FOR PRINCIPAL MSE                        #
#==================================================#
fn_principal_post<-function(DATA,POST,newtheta) {
  # INPUT   DATA (data.frame): data for one measure
  #         POST (list): posterior draws for one measure
  #         newtheta: the effect to be changed for future (* superscript)
  # OUTPUT  list: predicted future draws for eta and y
  R<-nrow(POST$eta)
  THETA<-POST$theta
  THETA[[1]]<-POST$theta[[1]]
  
  for (r in 2:length(POST$theta)) {
    for (s in 1:ncol(THETA[[r]])) {
      if (r %in% newtheta) {
        THETA[[r]][,s]<-rnorm(R,mean=0,sd=sqrt(POST$sigsq[[r]]))  
      } else {
        THETA[[r]]<-POST$theta[[r]]
      }
    }
  }
  
  # eta specification
  base<-vector('list')
  # effects
  base[[1]]<-'year'; base[[2]]<-'state'; base[[3]]<-'substate'
  # interaction effects
  # enter year first (sorting is based on year first)
  base[[4]]<-c('year','substate')
  
  param<-fn_idx(DATA=DATA,spec=base,theta0=F)
  
  ETA<-THETA[[1]][,param$idx[,1]]+
       THETA[[2]][,param$idx[,2]]+
       THETA[[3]][,param$idx[,3]]+
       THETA[[4]][,param$idx[,4]]

  ystar<-matrix(rnorm(nrow(DATA)*R,mean=t(ETA),sd=rep(sqrt(DATA$V),R)),nrow=nrow(DATA),ncol=R)
  
  return(list(theta=THETA,eta=ETA,sigsq=POST$sigsq,ystar=ystar))
}

#==================================================#
#   FUTURE DRAWS FOR ALL MEASURES                  #
#==================================================#
fn_get_pr_post<-function(measures,post) {
  # INPUT   measures (list): data frames for each measure
  #         post (list): posterior draws for each measure
  # OUTPUT  list: list of future draws for all measures
  pr_post<-list()
  for (m in names(measures)) {
    pr_post[[m]]<-fn_principal_post(DATA=measures[[m]],POST=post[[m]],newtheta=4)
  }
  
  return(pr_post)
}
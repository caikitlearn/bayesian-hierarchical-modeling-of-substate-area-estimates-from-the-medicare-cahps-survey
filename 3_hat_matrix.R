#==================================================#
#   ASSIGNING GROUPS TO HAT MATRIX                 #
#==================================================#
fn_Hgroups<-function(DATA) {
  # INPUT   DATA (data.frame): data for one measure
  # OUTPUT  matrix: assigned groups for each element of hat matrix
  Hgroups<-matrix(NA,nrow=nrow(DATA),ncol=nrow(DATA))
  for (i in 1:nrow(Hgroups)) {
    YEARi<-DATA$year[i]
    STATEi<-DATA$state[i]
    AREAi<-DATA$substate[i]
    for (j in 1:i) {
      YEARj<-DATA$year[j]
      STATEj<-DATA$state[j]
      AREAj<-DATA$substate[j]
      # same area, same year
      if (STATEi==STATEj & AREAi==AREAj & YEARi==YEARj) {
        Hgroups[i,j]<-1
        Hgroups[j,i]<-1
        # same area, different year
      } else if (STATEi==STATEj & AREAi==AREAj & YEARi!=YEARj) {
        Hgroups[i,j]<-2
        Hgroups[j,i]<-2
        # same state, different area, same year
      } else if (STATEi==STATEj & AREAi!=AREAj & YEARi==YEARj) {
        Hgroups[i,j]<-3
        Hgroups[j,i]<-3
        # same state, different area, different year
      } else if (STATEi==STATEj & AREAi!=AREAj & YEARi!=YEARj) {
        Hgroups[i,j]<-4
        Hgroups[j,i]<-4
        # different area, same year
      } else if (STATEi!=STATEj & AREAi!=AREAj & YEARi==YEARj) {
        Hgroups[i,j]<-5
        Hgroups[j,i]<-5
        # different area, different year
      } else if (STATEi!=STATEj & AREAi!=AREAj & YEARi!=YEARj) {
        Hgroups[i,j]<-6
        Hgroups[j,i]<-6
      }
    }
  }
  return(Hgroups)
}

#==================================================#
#   GENERATING HAT MATRIX DRAWS                    #
#==================================================#
fn_Hresults<-function(Hgroups,DATA,POST,param) {
  # INPUT   Hgroups (matrix): assigned groups for each element of hat matrix
  #         DATA (data.frame): data for one measure
  #         POST (list): posterior draws for one measure
  #         param (data.frame): indices from fn_idx()
  # OUTPUT  list: summaries of the hat matrix
  # calculating length of theta
  lentheta<-0
  for (i in 1:length(POST$theta)) {
    lentheta<-lentheta+ncol(POST$theta[[i]])
  }
  
  # constructing Y
  Y<-DATA$Ybar
  
  # constructing V
  V<-diag(DATA$V)
  
  # constructing X
  supps<-matrix(NA,nrow=nrow(param$idx),ncol=ncol(param$idx))
  supps[,1]<-param$idx[,1]
  for (i in 2:ncol(supps)) {
    add<-0
    for (j in 1:(i-1)) {
      add<-add+ncol(POST$theta[[j]])
    }
    supps[,i]<-param$idx[,i]+add
  }
  
  X<-matrix(0,nrow=length(Y),ncol=lentheta)
  for (i in 1:nrow(X)) {
    X[i,supps[i,]]<-1
  }
  
  # X has been checked so that X%*%theta.hat is equal to the overall posterior means
  R<-nrow(POST$eta)
  SIGMAS<-POST$sigsq
  THETAS<-POST$theta
  
  # posterior means
  MEANS<-matrix(NA,nrow=R,ncol=470)
  # results of the hat matrix
  hmean<-matrix(NA,nrow=R,ncol=6)
  hmed<-matrix(NA,nrow=R,ncol=6)
  hlow<-matrix(NA,nrow=R,ncol=6)
  hupp<-matrix(NA,nrow=R,ncol=6)
  
  for (t in 1:R) {
    if (t%%(R/100)==0) {cat(t/R*100,'% ',sep='')}
    
    s.diag<-NULL
    sigpar<-rep(NA,length(SIGMAS))
    
    for (s in 1:length(SIGMAS)) {
      sigpar[s]<-SIGMAS[[s]][t]
      s.diag<-c(s.diag,rep(sigpar[s],ncol(THETAS[[s]])))
    }
    
    OMEGA<-diag(s.diag)
    
    H<-X%*%solve(solve(OMEGA)+t(X)%*%solve(V)%*%X)%*%t(X)%*%solve(V)
    MEANS[t,]<-as.vector(H%*%DATA$Ybar)
    
    Hrowsums<-matrix(NA,nrow=nrow(H),ncol=6)
    for (r in 1:6) {
      for (i in 1:nrow(H)) {
        Hrowsums[i,r]<-sum(H[i,which(Hgroups[i,]==r)])
      }
    }
    
    hmean[t,]<-apply(Hrowsums,2,mean)
    hmed[t,]<-apply(Hrowsums,2,median)
    hlow[t,]<-apply(Hrowsums,2,quantile,probs=0.05)
    hupp[t,]<-apply(Hrowsums,2,quantile,probs=0.95)
  }
  
  return(list(MEANS=MEANS,hmean=hmean,hmed=hmed,hlow=hlow,hupp=hupp))
}

#==================================================#
#   SAVING HAT MATRIX DRAWS                        #
#==================================================#
fn_save_hat_matrix<-function() {
  measures<-fn_data_prep()
  
  # eta specification
  base<-vector('list')
  # effects
  base[[1]]<-'year'; base[[2]]<-'state'; base[[3]]<-'substate'
  # interaction effects
  # enter year first (sorting is based on year first)
  base[[4]]<-c('year','substate')
  
  # since all data have the same domains, all assigned groups for hat matrix are the same across measures
  Hgroups<-fn_Hgroups(DATA=measures[[names(measures)[1]]])
  
  for (m in names(measures)) {
    post<-readRDS(paste('posterior/base_post_',m,'.rds',sep=''))
    hat_matrix<-fn_Hresults(Hgroups=Hgroups,DATA=measures[[m]],POST=post,param=fn_idx(DATA=measures[[m]],spec=base,theta0=F))
    saveRDS(hat_matrix,paste('posterior/hat_',m,'.rds',sep=''))
  }
}
#==================================================#
#   INDEX GENERATION TO REPRESENT EFFECTS          #
#==================================================#
fn_idx<-function(DATA,spec,theta0=T) {
  # INPUT   DATA (data.frame): data for one measure
  #         spec (list): names of the fixed and random effects
  #         theta0 (boolean): if we use a theta0 fixed effect
  # OUTPUT  list: Tg: number of levels for each effect
  #               idx: indices matching data to each effect
  N<-nrow(DATA)
  
  # creating lengths of the random effects
  G<-length(spec)
  Tg<-rep(NA,G)
  
  # fixed and random effects
  for (g in 1:G) {
    # interaction effects
    if(length(spec[[g]])>1) {
      # find out lengths of the variables involved in the interaction
      int1<-length(table(subset(DATA,select=spec[[g]][1])))
      int2<-length(table(subset(DATA,select=spec[[g]][2])))
      # result is the product of the lengths
      Tg[g]<-int1*int2
    } else {
      Tg[g]<-length(table(subset(DATA,select=spec[[g]])))  
    }
  }
  
  # creating indices for each fixed and random effect
  idx<-matrix(nrow=N,ncol=G)
  
  # fixed and random effects
  for (g in 1:G) {
    # interaction effects
    if(length(spec[[g]])>1) {
      # find the variables involved in the interaction
      v1<-unique(subset(DATA,select=spec[[g]][1]))[,1]
      v2<-unique(subset(DATA,select=spec[[g]][2]))[,1]
      # create unique look-up terms for the interaction variables
      # all values of v2 are repeated for the first value of v1, and so on
      xterms<-paste(rep(v1,rep(length(v2),length(v1))),v2)
      idx[,g]<-match(paste(subset(DATA,select=spec[[g]][1])[,1],
                           subset(DATA,select=spec[[g]][2])[,1]),xterms)
    } else {
      # find the index for that variable
      idx[,g]<-match(subset(DATA,select=spec[[g]])[,1],
                     unique(subset(DATA,select=spec[[g]]))[,1])
    }
  }
  
  if (theta0) {
    Tg<-c(1,Tg)
    idx<-cbind(1,idx)
  }
  
  return(list(Tg=Tg,idx=idx))
}

#==================================================#
#   MCMC ALGORITHM                                 #
#==================================================#
fn_unimodel<-function(DATA,param,R=10000,burn=1000,exvar=0,tau=10^3,alpha=10^-3,beta=10^-3) {
  # INPUT   DATA (data.frame): data for one measure
  #         R (integer): number of iterations
  #         param (list): output from fn_idx
  #         burn (integer): burn-in period
  #         tau (number): variance of prior on fixed effect
  #         alpha (number): alpha of IG prior
  #         beta (number): beta of IG prior
  # OUTPUT  list: theta: posterior theta parameters
  #               sigsq: posterior sigma parameters
  Tg<-param$Tg
  idx<-param$idx
  
  R<-R+burn
  
  G<-length(Tg)
  
  if (length(exvar)==1) {
    exvar<-rep(0,G)
  }
  
  # create a list of matrices to store the posterior theta and sigma
  theta<-list()
  sigsq<-list()
  for (g in 1:G) {
    theta[[g]]<-matrix(NA,nrow=R,ncol=Tg[g])
    if (exvar[g]>0) {
      Sg<-Tg[exvar[g]]
    } else {
      Sg<-1
    }
    sigsq[[g]]<-matrix(NA,nrow=R,ncol=Sg)
  }
  
  eta<-matrix(NA,nrow=R,ncol=nrow(DATA))
  
  # initial values
  for (g in 1:G) {
    theta[[g]][1,]<-0
    sigsq[[g]][1,]<-0.05
  }
  sigsq[[1]][,1]<-tau^2 # fixed effect
  
  eta[1,]<-0
  
  yi<-DATA$Ybar
  vi<-DATA$V
  
  for (i in 2:R) {
    
    # progress bar
    if (i%%(R/100)==0) {cat(i/R*100,'% ',sep='')}
    
    # update theta_g
    for (g in 1:G) {
      # constructing the r_i
      ri<-yi
      for (j in (1:G)[-g]) {
        if (j<g) {
          ri<-ri-theta[[j]][i,idx[,j]]
        } else {
          ri<-ri-theta[[j]][i-1,idx[,j]]
        }
      }
      # drawing theta
      for (t in 1:Tg[g]) {
        Igt<-which(idx[,g]==t)
        # one sigsq for each effect
        if (ncol(sigsq[[g]])==1) {
          s_gs<-sigsq[[g]][i-1,]
        } else if (ncol(sigsq[[g]])>1) { # if sigsq varies by effect
          st<-idx[t,exvar[[g]]]
          s_gs<-sigsq[[g]][i-1,st]
        }
        var.theta<-(sum(1/vi[Igt])+1/s_gs)^-1  
        mean.theta<-sum((ri/vi)[Igt])*var.theta
        theta[[g]][i,t]<-rnorm(1,mean=mean.theta,sd=sqrt(var.theta))
      }
    }
    
    # update sigsq
    for (g in 2:G) {
      if (ncol(sigsq[[g]])==1) {
        shape.s<-alpha+Tg[g]/2
        rate.s<-beta+sum(theta[[g]][i,]^2)/2
        sigsq[[g]][i,]<-1/rgamma(1,shape=shape.s,rate=rate.s)
      } else if (ncol(sigsq[[g]])>1) {
        # multiple sigma parameters
        for (s in 1:Tg[exvar[g]]) {
          Tgs<-sort(unique(idx[which(idx[,exvar[g]]==s),g]))
          shape.s<-alpha+length(Tgs)/2
          rate.s<-beta+sum(theta[[g]][i,Tgs]^2)/2
          sigsq[[g]][i,s]<-1/rgamma(1,shape=shape.s,rate=rate.s)   
        } 
      }
    }
    
    # compute eta
    sums<-0
    for (g in 1:G) {
      sums<-sums+theta[[g]][i,idx[,g]]
    }
    eta[i,]<-sums
  }
  
  theta.b<-list()
  sigsq.b<-list()
  for (g in 1:G) {
    theta.b[[g]]<-theta[[g]][(burn+1):R,]
    sigsq.b[[g]]<-sigsq[[g]][(burn+1):R,]
  }
  
  return(list(theta=theta.b,
              eta=eta[(burn+1):R,],
              sigsq=sigsq.b))
}

#==================================================#
#   THINNING DATA                                  #
#==================================================#
fn_thin<-function(DATA,len=50) {
  # INPUT   DATA (vector, matrix, or array): chain to be thinned
  #         len (integer): amount of thinning
  # OUTPUT  vector, matrix, or array: thinned version of input
  
  # thinning vectors
  if (is.null(dim(DATA))) {
    thin<-seq(1,length(DATA),by=len)
    DATA<-DATA[thin]
    # thinning matrices
  } else if (length(dim(DATA))==2) {
    thin<-seq(1,nrow(DATA),by=len)
    DATA<-DATA[thin,]
    # thinning arrays of length 3
  } else if (length(dim(DATA))==3) {
    thin<-seq(1,dim(DATA)[3],by=len)
    DATA<-DATA[,,thin]
    # thinning arrays of length 4
  } else if (length(dim(DATA))==4) {
    thin<-seq(1,dim(DATA)[4],by=len)
    DATA<-DATA[,,,thin]
  }
  
  return(DATA)
}

#==================================================#
#   THINNING POSTERIOR                             #
#==================================================#
fn_post_thin<-function(POST,len=50) {
  # INPUT   POST (list): posterior parameters to be thinned
  #         len (integer): amount of thinning
  # OUTPUT  list: thinned version of input
  theta<-list()
  sigsq<-list()
  
  for (i in 1:length(POST$theta)) {
    theta[[i]]<-fn_thin(DATA=POST$theta[[i]],len=len)
  }
  
  for (i in 1:length(POST$sigsq)) {
    sigsq[[i]]<-fn_thin(DATA=POST$sigsq[[i]],len=len)
  }
  
  return(list(theta=theta,eta=fn_thin(DATA=POST$eta,len=len),sigsq=sigsq))
}

#==================================================#
#   SAVING SIMULATIONS                             #
#==================================================#
fn_save_post<-function(R=50000,len=50) {
  # INPUT   R (int): number of posterior draws
  #         len (int): amount of thinning
  measures<-fn_data_prep()
  
  # eta specification
  # for interaction effects, enter year first (sorting is based on year first)
  annu<-list()
  annu[[1]]<-'state';annu[[2]]<-'substate'
  
  base<-list()
  base[[1]]<-'year';base[[2]]<-'state';base[[3]]<-'substate';base[[4]]<-c('year','substate')
  
  sby<-list()
  sby[[1]]<-'year';sby[[2]]<-'state';sby[[3]]<-'substate';sby[[4]]<-c('year','state');sby[[5]]<-c('year','substate')
  
  rich<-list()
  rich[[1]]<-'year';rich[[2]]<-'state';rich[[3]]<-'substate';rich[[4]]<-c('year','state');rich[[5]]<-c('year','substate')
  
  full<-list()
  full[[1]]<-'year';full[[2]]<-'state';full[[3]]<-'substate';full[[4]]<-c('year','state');full[[5]]<-c('year','substate')
  
  # simulating the posterior distributions
  for (out in names(measures)) {
    cat('measure: ',out,' ',sep='')
    
    cat('\nmodel: annu ')
    y_terms<-list()
    years<-2012:2016
    for (y in 1:length(years)) {
      cat('\nyear: ')
      cat(years[y])
      cat(' ')
      y.data<-measures[[out]][measures[[out]]$year==years[y],]
      param.annu<-fn_idx(DATA=y.data,spec=annu)
      m<-fn_unimodel(DATA=y.data,param=param.annu,R=R)
      y_terms[[y]]<-fn_post_thin(POST=m,len=len)
    }
    annu_post<-y_terms
    saveRDS(annu_post,paste('posterior/annu_post_',out,'.rds',sep=''))
    
    cat('\nmodel: base ')
    param.base<-fn_idx(DATA=measures[[out]],spec=base,theta0=F)
    base_post<-fn_post_thin(POST=fn_unimodel(DATA=measures[[out]],param=param.base,R=R),len=len)
    saveRDS(base_post,paste('posterior/base_post_',out,'.rds',sep=''))
    
    cat('\nmodel: sby ')
    param.sby<-fn_idx(DATA=measures[[out]],spec=sby,theta0=F)
    sby_post<-fn_post_thin(POST=fn_unimodel(DATA=measures[[out]],param=param.sby,R=R),len=len)
    saveRDS(sby_post,paste('posterior/sby_post_',out,'.rds',sep=''))
    
    cat('\nmodel: rich ')
    param.rich<-fn_idx(DATA=measures[[out]],spec=rich,theta0=F)
    rich_post<-fn_post_thin(POST=fn_unimodel(DATA=measures[[out]],param=param.rich,R=R,exvar=c(0,0,2,0,0)),len=len)
    saveRDS(rich_post,paste('posterior/rich_post_',out,'.rds',sep=''))
    
    cat('\nmodel: full ')
    param.full<-fn_idx(DATA=measures[[out]],spec=full,theta0=T)
    full_post<-fn_post_thin(POST=fn_unimodel(DATA=measures[[out]],param=param.full,R=R),len=len)
    saveRDS(full_post,paste('posterior/full_post_',out,'.rds',sep=''))
    
    cat('\n')
  }
}
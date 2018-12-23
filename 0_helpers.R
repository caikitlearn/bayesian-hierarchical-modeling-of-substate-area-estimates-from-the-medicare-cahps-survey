#==================================================#
#   LOADING DATA                                   #
#==================================================#
fn_data_prep<-function(filename='data/raw_ffs.csv') {
  # INPUT   filename (string): file name for data csv file
  # OUTPUT  list: data frames for each measure
  # full data
  ffs_full<-read.table(filename,sep=',',header=T)
  ffs_full$Ybar<-ffs_full$score
  ffs_full$V<-ffs_full$vp
  
  # extracting variables for analysis
  ffs<-subset(ffs_full,select=c(item,state,substate,year,Ybar,V,usen))

  measures<-list()
  for (m in levels(ffs$item)) {
    measures[[m]]<-ffs[ffs$item==m,]
  }
  
  return(measures)
}

#==================================================#
#   LOADING SAVED POSTERIOR DRAWS                  #
#==================================================#
fn_load_post<-function(measures) {
  # INPUT   measures (list): data frames for each measure
  # OUTPUT  list: posterior draws for each measure
  post<-list()
  for (m in names(measures)) {
    post[[m]]<-readRDS(paste('posterior/base_post_',m,'.rds',sep=''))
  }
  return(post)
}

#==================================================#
#   PRINTS RESULTS FOR TABLE 3                     #
#==================================================#
fn_table3<-function(model_comp) {
  # INPUT   model_comp (list): results from model comparison
  print(round(model_comp$CVerrs,1))
  print(t(apply(model_comp$CVerrs,1,function(xx){rank(-xx)})))
}

#==================================================#
#   PRINTS RESULTS FOR TABLE 4                     #
#==================================================#
fn_table4<-function(model_comp) {
  # INPUT   model_comp (list): results from model comparison
  print(model_comp$PV)
}

#==================================================#
#   PRINTS RESULTS FOR TABLE 5                     #
#==================================================#
fn_table5<-function() {
  hat_gcq_comp<-readRDS('posterior/hat_gcq_comp.rds')
  hat_cs_comp<-readRDS('posterior/hat_cs_comp.rds')
  t5_gcq<-t(rbind(round(colMeans(hat_gcq_comp$hmean),2),round(apply(hat_gcq_comp$hmean,2,quantile,probs=c(0.05,0.95)),2)))
  t5_cs<-t(rbind(round(colMeans(hat_cs_comp$hmean),2),round(apply(hat_cs_comp$hmean,2,quantile,probs=c(0.05,0.95)),2)))
  print(cbind(t5_gcq,t5_cs))
}

#==================================================#
#   STANDARDIZING VECTORS                          #
#==================================================#
fn_std<-function(xx) {return((xx-mean(xx))/sd(xx))}
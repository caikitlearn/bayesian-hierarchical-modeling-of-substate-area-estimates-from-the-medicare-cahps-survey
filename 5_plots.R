#==================================================#
#   QQ PLOT FOR A PARTICULAR EFFECT                #
#==================================================#
fn_plot_qq<-function(POST,N,title) {
  # INPUT   POST (list): posterior draws for one measure
  #         N (int): number of draws to sample
  #         title: subplot title
  std_thinned<-as.vector(apply(POST,2,fn_std))
  std_thinned<-sample(std_thinned,size=N)
  qqnorm(std_thinned,xlim=c(-4,4),ylim=c(-6,6),main=title,pch='.',cex.main=2,cex.axis=2)
  
  abline(0,1,col='red',lty=2)  
}

#==================================================#
#   FIGURE 1: QQ PLOTS FOR RANDOM EFFECTS          #
#==================================================#
fn_figure1<-function(N=5000) {
  # INPUT   N (int): number of draws to sample
  bitmap('figures/fig1.tiff',height=4,width=9,units='in',type='tifflzw',res=800)
  post_rate_plan<-readRDS('posterior/base_post_rate_plan.rds')  
  
  lay<-matrix(c(1,1,1,2,3,4),nrow=2,ncol=3,byrow=T)
  layout(mat=lay,heights=c(0.1,0.9))
  par(mar = c(0,0,3,0))
  plot(1,type='n',axes=F,xlab='',ylab='')
  
  par(mar = c(5.1,4.1,4.1,2.1))
  fn_plot_qq(POST=post_rate_plan$theta[[2]],N=N,title='State Random Effects')
  fn_plot_qq(POST=post_rate_plan$theta[[3]],N=N,title='Substate Random Effects')
  fn_plot_qq(POST=post_rate_plan$theta[[4]],N=N,title='Substate by Year Random Effects')
  
  dev.off()
}

#==================================================#
#   SUMMARIES FOR FIGURE 2                         #
#==================================================#
fn_summary_2<-function(DATA,POST) {
  # INPUT   DATA (data.frame): data for one measure
  #         POST (list): posterior draws for one measure
  # OUTPUT  data.frame: mean, median, and 25/75 quantiles
  OUT<-matrix(NA,nrow=4,ncol=length(POST$sigsq)-1)
  postmean<-rep(NA,length(POST$sigsq)-1)
  post25<-rep(NA,length(POST$sigsq)-1)
  post50<-rep(NA,length(POST$sigsq)-1)
  post75<-rep(NA,length(POST$sigsq)-1)
  for (i in 2:length(POST$sigsq)) {
    postmean[i-1]<-round(mean(POST$sigsq[[i]])/mean(DATA$V),2)
    post25[i-1]<-round(quantile(POST$sigsq[[i]]/mean(DATA$V),probs=c(0.25)),2)
    post50[i-1]<-round(quantile(POST$sigsq[[i]]/mean(DATA$V),probs=c(0.5)),2)
    post75[i-1]<-round(quantile(POST$sigsq[[i]]/mean(DATA$V),probs=c(0.75)),2)
  }
  
  OUT[1,]<-post25
  OUT[2,]<-post50
  OUT[3,]<-postmean
  OUT[4,]<-post75
  return(OUT)
}

#==================================================#
#   FIGURE 2: RANDOM EFFECT VARIANCES              #
#==================================================#
fn_figure2<-function(measures,post) {
  # INPUT   measures (list): data frames for each measure
  #         post (list): posterior draws for each measure
  # OUTPUT  saves figure2 to /figures
  lbls<-c('Coordination','of Care','Customer','Service','Get Care','Quickly',
          'Get Need-','ed Care','Flu Immu-','ization',
          'Doctor Com-','munication','Rating','of Care','Rating','of Plan')
  
  xloc<-sort(c(seq(0.8,7.8,1),seq(1,8,1),seq(1.2,8.2,1)),decreasing=T)
  Ypt<-NULL
  for (m in names(measures)) {
    Ypt<-cbind(Ypt,fn_summary_2(DATA=measures[[m]],POST=post[[m]]))
  }
  
  clrs<-'black'
  
  lvs<-rep(c('S','A','AY'),8)
  
  xmin<-0.045
  xmax<-11
  
  bitmap('figures/fig2.tiff',height=5,width=5,units='in',type='tifflzw',res=800)
  
  par(mar=c(6.1,6.1,6.1,1.1),xpd=T)
  plot(Ypt[4,],xloc,log='x',col=clrs,xlab=expression(hat(sigma)^2/bar(v)),
       yaxt='n',ylab='',pch='|',xlim=c(xmin,xmax))
  mtext('S: State; A: Substate Area; AY: Substate Area by Year',side=3,
        at=c(0.6),line=0.7)
  legend('top',inset=c(0,-0.12),legend=c('Posterior Mean','Posterior Median','25-75%ile'),pch=c(2,1,45),
         bty='n',horiz=T)
  text(Ypt[4,]*1.2,xloc,labels=lvs,col=clrs)
  axis(2,at=c(8.2,7.8,
              7.2,6.8,
              6.2,5.8,
              5.2,4.8,
              4.2,3.8,
              3.2,2.8,
              2.2,1.8,
              1.2,0.8),labels=lbls,las=2,cex.axis=1,tck=0)
  points(Ypt[3,],xloc,col=clrs,pch=2)
  points(Ypt[2,],xloc,col=clrs,cex=1)
  points(Ypt[1,],xloc,col=clrs,pch='|',cex=1)
  for (i in 1:length(xloc)) {
    lines(c(Ypt[1,i],Ypt[4,i]),rep(xloc[i],2),lty=1,lwd=0.5)  
  }
  
  linesy<-seq(1.5,7.5,by=1)
  for (i in 1:length(linesy)) {
    lines(c(xmin,xmax),rep(linesy[i],2),col='grey',lty=3)  
  }
  
  dev.off()  
}

#==================================================#
#   SUMMARIES FOR FIGURE 3                         #
#==================================================#
fn_summary_3<-function(DATA,POST) {
  # INPUT   DATA (data.frame): data for one measure
  #         POST (list): posterior draws for one measure
  # OUTPUT  data.frame: mean, median, and 25/75 quantiles
  OUT<-matrix(NA,nrow=4,ncol=3)
  eta.var<-apply(POST$eta,2,var)
  vr<-eta.var/DATA$V
  
  areas<-as.vector(table(DATA$state)/dim(table(DATA$year)))
  states<-names(table(DATA$state))
  
  # group by states with 2, 3-4, 5-7 substates
  areasizes<-areas[match(DATA$state,states)]
  smalls<-vr[areasizes==2]
  mids<-vr[areasizes==3|areasizes==4]
  larges<-vr[areasizes>4]
  
  OUT[1,]<-c(quantile(smalls,0.25),quantile(mids,0.25),quantile(larges,0.25))
  OUT[2,]<-c(median(smalls),median(mids),median(larges))
  OUT[3,]<-c(mean(smalls),mean(mids),mean(larges))
  OUT[4,]<-c(quantile(smalls,0.75),quantile(mids,0.75),quantile(larges,0.75))
  return(OUT)
}

#==================================================#
#   FIGURE 3: VARIANCE REDUCTIONS                  #
#==================================================#
fn_figure3<-function(measures,post) {
  # INPUT   measures (list): data frames for each measure
  #         post (list): posterior draws for each measure
  # OUTPUT  saves figure3 to /figures
  lbls<-c('Coordination','of Care','Customer','Service','Get Care','Quickly',
          'Get Need-','ed Care','Flu Immu-','ization',
          'Doctor Com-','munication','Rating','of Care','Rating','of Plan')
  
  Ypt<-NULL
  for (m in names(measures)) {
    Ypt<-cbind(Ypt,fn_summary_3(DATA=measures[[m]],POST=post[[m]]))
  }
  
  clrs<-'black'
  lvs<-rep(c('S','M','L'),8)
  
  xloc<-c(0.8,1,1.2)+rep(seq(0,7,1),rep(3,8))
  xloc<-rev(xloc)
  
  bitmap('figures/fig3.tiff',height=5,width=5,units='in',type='tifflzw',res=800)
  
  par(mar=c(6.1,6.1,6.1,1.1),xpd=T)
  
  xmin<-0.1
  xmax<-0.6
  
  plot(Ypt[4,],xloc,col=clrs,xlab='Posterior Variance / Sampling Variance',
       yaxt='n',ylab='',pch='|',xlim=c(xmin,xmax))
  mtext('S: 2 substates; M: 3-4 substates; L: 5+ substates',side=3,
        at=c(0.335),line=0.7)
  legend('top',inset=c(0,-0.12),legend=c('Posterior Mean','Posterior Median','25-75%ile'),pch=c(2,1,45),
         bty='n',horiz=T)
  text(Ypt[4,]+0.02,xloc,labels=lvs,col=clrs)
  axis(2,at=c(8.2,7.8,
              7.2,6.8,
              6.2,5.8,
              5.2,4.8,
              4.2,3.8,
              3.2,2.8,
              2.2,1.8,
              1.2,0.8),labels=lbls,las=2,cex.axis=1,tck=0)
  points(Ypt[3,],xloc,col=clrs,pch=2) # post mean
  points(Ypt[2,],xloc,col=clrs,cex=1)
  points(Ypt[1,],xloc,col=clrs,pch='|',cex=1)
  for (i in 1:length(xloc)) {
    lines(c(Ypt[1,i],Ypt[4,i]),rep(xloc[i],2),lty=1,lwd=0.5)  
  }
  
  linesy<-seq(1.5,7.5,by=1)
  for (i in 1:length(linesy)) {
    lines(c(xmin,xmax),rep(linesy[i],2),col='grey',lty=3)  
  }
  
  dev.off()
}

#==================================================#
#   MSE CALCULATION                                #
#==================================================#
fn_calc_mse<-function(DATA,POST) {
  # INPUT   DATA (data.frame): data for one measure
  #         POST (list): posterior draws for one measure
  # OUTPUT  vector: mse ratios for each state
  states<-names(table(DATA$state))
  years<-as.numeric(names(table(DATA$year)))
  
  eta<-POST$eta
  ystar<-POST$ystar
  
  postmean<-apply(eta,2,mean)
  postvar<-apply(eta,2,var)
  
  MSE.un<-matrix(NA,ncol=length(states),nrow=ncol(ystar))
  MSE.po<-matrix(NA,ncol=length(states),nrow=ncol(ystar))
  
  for (r in 1:ncol(ystar)) {
    terms<-as.data.frame(cbind(DATA,postmean,postvar,eta[r,]))
    for (i in 1:length(states)) {
      subd<-terms[terms$year==years[5]&terms$state==states[i],]
      MSE.un[r,i]<-sum(subd$V)
      statevar<-sum(subd$V*subd$usen^2)/(sum(subd$usen)^2)
      etabar<-sum((subd$`eta[r, ]`)*subd$usen)/sum(subd$usen)
      MSE.po[r,i]<-nrow(subd)*statevar+sum((etabar-subd$`eta[r, ]`)^2)
    }
  }
  
  return(colMeans(MSE.po)/(colMeans(MSE.po)+colMeans(MSE.un)))
}

#==================================================#
#   PLOTTING MSE FOR ONE MEASURE                   #
#==================================================#
fn_plot_mse<-function(MSEratio,states,areas,nn) {
  # INPUT   MSEratio (vector): mse ratios for each state
  #         states (vector): vector of string abbreviations for each state
  #         areas (vector): vector of number of areas per state
  #         nn (int): index for plotting title
  ordered<-cbind(areas,MSEratio)
  MSEsorted<-ordered[order(ordered[,1]),]

  n.unpooled<-sum(MSEratio>0.5)
  n.pooled<-sum(MSEratio<=0.5)

  NAMES<-c('Coordination of Care','Customer Service','Get Care Quickly','Get Needed Care','Flu Immunization',
           'Doctor Communication','Rating of Care','Rating of Plan')

  ymin<-0
  ymax<-1

  I<-length(MSEratio)

  plot(1:I,MSEsorted[,2],xaxt='n',yaxt='n',ylab='MSE Ratio',
       ylim=c(ymin,ymax),xlab='',main=NAMES[nn],cex.lab=1.25,cex.main=1.5)
  
  lines(c(0.5,32.5),c(0.5,0.5),col='grey')
  axis(1,at=1:I,labels=rownames(MSEsorted),las=2,cex.axis=1)
  axis(2,at=c(0.2,0.4,0.5,0.6,0.8))
  text(30,0.92,n.unpooled,cex=1.25)
  text(10,0.92,paste('Substate-level Reporting'),cex=1.25)
  text(30,0.08,n.pooled,cex=1.25)
  text(10,0.08,paste('State-level Reporting'),cex=1.25)
}


#==================================================#
#   FIGURE 4: MSE PLOTS                            #
#==================================================#
fn_figure4<-function(measures,post) {
  # INPUT   measures (list): data frames for each measure
  #         post (list): posterior draws for each measure
  pr_post<-fn_get_pr_post(measures=measures,post=post)
  
  mse<-list()
  for (m in names(measures)) {
    mse[[m]]<-fn_calc_mse(DATA=measures[[m]],POST=pr_post[[m]])
  }
  
  bitmap('figures/fig4.tiff',height=7,width=5,units='in',type='tifflzw',res=800)
  
  lay<-matrix(c(1,1,2,3,4,5,6,7,8,9),nrow=5,ncol=2,byrow=T)
  layout(mat=lay,heights=c(0.06,0.235,0.235,0.235,0.235))
  par(mar = c(0,0,3,0))
  plot(1,type='n',axes=F,xlab='',ylab='')
  
  par(mar = c(3.1,3.1,3.1,3.1))
  states<-levels(measures[[names(measures)[[1]]]]$state)
  areas<-table(measures[[names(measures)[[1]]]]$state)/5
  for (i in 1:length(mse)) {
    fn_plot_mse(MSEratio=mse[[names(mse)[i]]],states=states,areas=areas,nn=i)
  }
  
  dev.off()  
}
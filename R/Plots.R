firingRateMapPlot <- function(m,name="",
                                 outma=c(2.0,2.0,2.0,2.0),margin=c(1,1,1,1),
                                 axis.x.mgp=c(0,0,0),axis.y.mgp=c(0,0,0),
                                 cex.x.axis=0.6,cex.y.axis=0.6,cex.lab=0.6,
                                 xlab="",ylab="",show.xlab=TRUE,main.title="",peak.rate.prefix="")
{
  jet.colors = colorRampPalette(c("#00007F", "blue","#007FFF",  "cyan", "#7FFF7F", "yellow", "#FF7F00","red"))
  par(oma=outma,mar=margin)
  
  if(class(m)=="matrix"){
    image(m,zlim=c(0,max(m,na.rm=T)), col=jet.colors(200),xlab='',ylab='',axes=FALSE)
    mtext(paste(peak.rate.prefix,name,round(max(m,na.rm=T),digits=2),"Hz"),line=-0.1,cex=0.6,side=3)
    
  }
  if(class(m)=="data.frame"){
    df<-m
    xlen <- length(unique(df$x))
    ylen <- length(unique(df$y))
    zzz <- matrix(df$rate,nrow=ylen,ncol=xlen)
    image(unique(df$y),unique(df$x),zzz,zlim=c(0,max(df$rate,na.rm=T)), col=jet.colors(200),xlab='',ylab='',axes=FALSE)
    mtext(paste(peak.rate.prefix,round(max(df$rate,na.rm=T),digits=2),"Hz"),side=3,at=median(unique(df$x)),line=-0.1,cex=0.6)
  }
  if(main.title!="")
  {
    mtext(main.title,side=3,line=0.3,cex=0.5)
  }
}


firingRateMapsPlot<-function(maps,names,fn="page.full.plot.pdf"){
  num.cols<-5
  num.rows<-6
  plot.per.page=num.cols*num.rows
  m<-matrix(c(rep(seq(0,1-(1/num.cols),1/num.cols),num.rows),
              rep(seq(1/num.cols,1,1/num.cols),num.rows),
              rep(seq(1-(1/num.rows),0,0-1/num.rows),each=num.cols),
              rep(seq(1,1/num.rows,0-1/num.rows),each=num.cols)),ncol=4)
#  pdf(file=fn,onefile=TRUE,paper="a4",width=8,height=10)
  index=1
  nCells<-dim(maps)[3]
  for (i in 1:nCells){
    if(index==1)
    {
      split.screen(m)  
    }
    screen(index)
    ## insert your plot function here
    firingRateMapPlot(maps[,,i],name=names[i])
    if(index==plot.per.page)
    {
      close.screen( all.screens = TRUE )
      index=0
    }
    index=index+1
  }
  close.screen(all.screens = TRUE)
#  dev.off()
}
firingRateMapAutoPlot <- function(m,name="",
                              outma=c(2.0,2.0,2.0,2.0),margin=c(1,1,1,1),
                              axis.x.mgp=c(0,0,0),axis.y.mgp=c(0,0,0),
                              cex.x.axis=0.6,cex.y.axis=0.6,cex.lab=0.6,
                              xlab="",ylab="",show.xlab=TRUE,main.title="",peak.rate.prefix="")
{
  jet.colors = colorRampPalette(c("#00007F", "blue","#007FFF",  "cyan", "#7FFF7F", "yellow", "#FF7F00","red"))
  par(oma=outma,mar=margin)
  image(m,zlim=c(-1,1), col=jet.colors(200),xlab='',ylab='',axes=FALSE)
  mtext(paste(peak.rate.prefix,name,round(max(m,na.rm=T),digits=2),"Hz"),line=-0.1,cex=0.6,side=3)
  if(main.title!="")
  {
    mtext(main.title,side=3,line=0.3,cex=0.5)
  }
}

firingRateMapAutosPlot<-function(maps,names,fn="page.full.plot.pdf"){
  num.cols<-5
  num.rows<-6
  plot.per.page=num.cols*num.rows
  m<-matrix(c(rep(seq(0,1-(1/num.cols),1/num.cols),num.rows),
              rep(seq(1/num.cols,1,1/num.cols),num.rows),
              rep(seq(1-(1/num.rows),0,0-1/num.rows),each=num.cols),
              rep(seq(1,1/num.rows,0-1/num.rows),each=num.cols)),ncol=4)
  #  pdf(file=fn,onefile=TRUE,paper="a4",width=8,height=10)
  index=1
  nCells<-dim(maps)[3]
  for (i in 1:nCells){
    if(index==1)
    {
      split.screen(m)  
    }
    screen(index)
    ## insert your plot function here
    firingRateMapAutoPlot(maps[,,i],name=names[i])
    if(index==plot.per.page)
    {
      close.screen( all.screens = TRUE )
      index=0
    }
    index=index+1
  }
  close.screen(all.screens = TRUE)
  #  dev.off()
}


linearRatePlot<-function(sp1,n=1,axis.y.pos=0,axis.x.pos=0,axis.y.las=2,
                           mgp.x=c(0.5,0.05,0.0),mgp.y=c(.8,0.3,0.2),xlab="Position (cm)",ylab="Rate (Hz)",
                           outma=c(1,1,1,0),margin=c(1.5,2,0.8,0.3),
                           xaxis.at=seq(0,80,20),yaxis.at=seq(0,50,10),...)
{
  rate=sp1@rateHisto[,n]
  rate[which(rate==-1.0)]<-NA
  position = seq(sp1@cmPerBin,sp1@cmPerBin*sp1@nBinRateHisto,sp1@cmPerBin)
  plotxlim=c(0,max(position))
  plotylim=c(0,max(rate,na.rm=T))
  par(mar=margin, oma=outma,cex.lab=0.6,cex.axis=0.6)
  plot (x=plotxlim, y=plotylim,type='n', axes=FALSE, pch=20,lwd=1,xlab="",ylab="")
  axis(side = 1, pos=axis.x.pos, at=xaxis.at, tck=-0.05,cex.axis=0.60,mgp=mgp.x)
  par(mgp=mgp.y)
  axis(side = 2, las=axis.y.las, pos=axis.y.pos,tck=-0.05,cex.axis=0.60,mgp=mgp.y)
  lines (position,rate,type='l',pch=20,xlab='',ylab='',lwd=0.75,col="black")
  title(xlab=xlab,mgp=mgp.x)
  title(ylab=ylab,mgp=mgp.y)
}


linearRatePlots<-function(sp1,fn="page.full.plot.pdf"){
  num.cols<-5
  num.rows<-6
  plot.per.page=num.cols*num.rows
  m<-matrix(c(rep(seq(0,1-(1/num.cols),1/num.cols),num.rows),
              rep(seq(1/num.cols,1,1/num.cols),num.rows),
              rep(seq(1-(1/num.rows),0,0-1/num.rows),each=num.cols),
              rep(seq(1,1/num.rows,0-1/num.rows),each=num.cols)),ncol=4)
  #  pdf(file=fn,onefile=TRUE,paper="a4",width=8,height=10)
  index=1
  for (i in 1:length(sp1@cellList)){
    if(index==1)
    {
      split.screen(m)  
    }
    screen(index)
    ## insert your plot function here
    linearRatePlot(sp1,i)
    if(index==plot.per.page)
    {
      close.screen( all.screens = TRUE )
      index=0
    }
    index=index+1
  }
  close.screen(all.screens = TRUE)
  #  dev.off()
}


firingRateMapPlot <- function(m,name="",
                                 outma=c(2.0,2.0,2.0,2.0),margin=c(1,1,1,1),
                                 axis.x.mgp=c(0,0,0),axis.y.mgp=c(0,0,0),
                                 cex.x.axis=0.6,cex.y.axis=0.6,cex.lab=0.6,
                                 xlab="",ylab="",show.xlab=TRUE,main.title="",peak.rate.prefix="")
{
  jet.colors = colorRampPalette(c("#00007F", "blue","#007FFF",  "cyan", "#7FFF7F", "yellow", "#FF7F00","red"))
  par(oma=outma,mar=margin)
  image(m,zlim=c(0,max(m,na.rm=T)), col=jet.colors(200),xlab='',ylab='',axes=FALSE)
  mtext(paste(peak.rate.prefix,name,round(max(m,na.rm=T),digits=2),"Hz"),line=-0.1,cex=0.6,side=3)
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
      close.screen( all = TRUE )
      index=0
    }
    index=index+1
  }
  close.screen(all = TRUE)
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
      close.screen( all = TRUE )
      index=0
    }
    index=index+1
  }
  close.screen(all = TRUE)
  #  dev.off()
}

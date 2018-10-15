#' Plot a spikes on the animal path
#' 
#' Plot the animal path from a Positrack object and 
#' the spikes in space from a SpatialProperties2d object.
#' 
#' 
#' @param sop List containing the data to plot the spike on path, you can get this list with the function spikeOnPath
#' It contains 4 objects: xSpike, ySpike, cluSpike, xPath and yPath 
#' @param clu Clu number of the cell of interest.
#' @param name String with the name of the cell
#' @param outma Outer margins of the figure
#' @param margin Inner margins of the figure
#' @param cex.name Size of the font use for the name of the map
#' @param cex.point Size of the spike points
#' @param plotxlim Limits of the x axis
#' @param plotylim Limits of the y axis
#' @param mgp.x mgp for the x axis
#' @param mgp.y mgp for the y axis
#' @param axis.x.pos Y position of the x axis
#' @param axis.y.pos X position of the y axis
#' @param plot.axis Whether to plot the axes or not (TRUE or FALSE)
#' @param xlab Name to display under the x axis
#' @param ylab Name to display at the left of the y axis
#' @param col The color of the spikes
spikeOnPathPlot <- function(sop,clu,name="",
                        outma=c(2.0,2.0,2.0,2.0),margin=c(1.5,1.5,1,1),
                        cex.name=0.6,cex.point=0.5,
                        plotxlim=c(0,80),plotylim=c(0,80),
                        mgp.x=c(0.5,0.05,0.0),mgp.y=c(.8,0.3,0.2),
                        axis.x.pos=0,axis.y.pos=0,
                        plot.axis=TRUE,
                        xlab="",ylab="",
                        col="red")
{
  par(mar=margin, oma=outma,cex.lab=0.6,cex.axis=0.6)
  plot (x=plotxlim, y=plotylim,type='n', axes=FALSE, pch=20,lwd=1,xlab="",ylab="")
  lines(sop$xPath,sop$yPath,type='l',xlab=xlab,ylab=ylab)
  points(sop$xSpike[which(sop$cluSpike==clu)],sop$ySpike[which(sop$cluSpike==clu)],col=col,
         pch=20,cex=cex.point)
  if(plot.axis){
    axis(side = 1, pos=axis.x.pos,  tck=-0.05,cex.axis=0.60,mgp=mgp.x)
    axis(side = 2, las=2, pos=axis.y.pos, tck=-0.05,cex.axis=0.60,mgp=mgp.y)
  }
  title(xlab=xlab,mgp=mgp.x)
  title(ylab=ylab,mgp=mgp.y)
  if(name!="")
  {
    mtext(name,side=3,line=0.3,cex=cex.name)
  }
}


#' Plot the animal path
#' 
#' Plot the animal path from a Positrack object
#' 
#' @param pt Positrack object
#' @param name String with the name of the plot
#' @param outma Outer margins of the figure
#' @param margin Inner margins of the figure
#' @param cex.name Size of the font use for the name of the map
#' @param cex.line Size of the line
#' @param cex.lab Size of labels
#' @param cex.axis Size of axis
#' @param plotxlim Limits of the x axis
#' @param plotylim Limits of the y axis
#' @param mgp.x mgp for the x axis
#' @param mgp.y mgp for the y axis
#' @param axis.x.pos Y position of the x axis
#' @param axis.y.pos X position of the y axis
#' @param plot.axis Whether to plot the axes or not (TRUE or FALSE)
#' @param xlab Name to display under the x axis
#' @param ylab Name to display at the left of the y axis
#' @param col Color of the path
pathPlot <- function(pt,name="",
                     outma=c(2.0,2.0,2.0,2.0),margin=c(1.5,1.5,1,1),
                     cex.name=0.6,cex.line=0.2,cex.lab=0.6,cex.axis=0.6,
                     plotxlim=c(0,80),plotylim=c(0,80),
                     mgp.x=c(0.5,0.05,0.0),mgp.y=c(.8,0.3,0.2),
                     axis.x.pos=0,axis.y.pos=0,
                     plot.axis=TRUE,
                     xlab="",ylab="",
                     col="black")
{
  if(class(pt)!="Positrack")
    stop("pt is not of the Positrack class")
  if(length(pt@x)==0)
    stop("length(pt@x)==0")
  par(mar=margin, oma=outma,cex.lab=cex.lab,cex.axis=cex.axis)
  plot (x=plotxlim, y=plotylim,type='n', axes=FALSE, pch=20,lwd=1,xlab="",ylab="",cex=cex.line)
  lines(pt@x,pt@y,type='l',xlab=xlab,ylab=ylab,lwd=0.2,col=col)
  if(plot.axis){
    axis(side = 1, pos=axis.x.pos,  tck=-0.05,cex.axis=cex.axis,mgp=mgp.x)
    axis(side = 2, las=2, pos=axis.y.pos, tck=-0.05,cex.axis=cex.axis,mgp=mgp.y)
  }
  title(xlab=xlab,mgp=mgp.x)
  title(ylab=ylab,mgp=mgp.y)
  if(name!="")
  {
    mtext(name,side=3,line=0.3,cex=cex.name)
  }
}








#' Plot a color scale for the firing rate maps
#' 
#' Plot a color scale that can be used beside the plots obtained with firingRateMapPlot()
#' 
#' @param outma Outer margins of the figure
#' @param margin Inner margins of the figure
#' @param cex.labels Size of the font use for 0 and Max
colorScalePlot<-function(outma=c(0.2,0.2,0.2,0.2),margin=c(0,1,0,1),cex.labels=0.5){
  par(oma=outma,mar=margin)
  jet.colors = grDevices::colorRampPalette(c("#00007F", "blue","#007FFF",  "cyan", "#7FFF7F", "yellow", "#FF7F00","red"))
  m<-matrix(0:200,ncol = 1)
  graphics::image(m,zlim=c(0,max(m,na.rm=T)), col=jet.colors(200),xlab='',ylab='',axes=F)
  graphics::axis(side = 1,labels = F,tck=0.2)
  graphics::axis(side = 3,labels = F,tck=0.2)
  graphics::mtext("0",side=2,las=1,cex=cex.labels,line=0.2)
  graphics::mtext("Max",side=4,las=1,cex=cex.labels,line=0.2)
}


#' Plot a single firing rate map
#' 
#' Plot a 2-dimensional representation of firing rate using the image function included in the graphics package
#' A square red map is printed if all values in the map are invalid
#' 
#' @param m A data.frame or matrix containing the data of the firing rate map. If a data.frame is given as argument, it should have the column names x, y and rate
#' @param name Character vector giving the name of the firing rate map which is displayed before the max firing rate
#' @param outma Outer margins of the figure
#' @param margin Inner margins of the figure
#' @param cex.title Size of the font use for the title
#' @param cex.name Size of the font use for the name of the map
#' @param xlab Name to display under the x axis
#' @param ylab Name to display at the left of the y axis
#' @param main.title A title for the figure
#' @param peak.rate.prefix Additional information to display before the peak firing rate.
firingRateMapPlot <- function(m,name="",
                              outma=c(2.0,2.0,2.0,2.0),margin=c(1,1,1,1),
                              cex.title=0.5,cex.name=0.5,
                              xlab="",ylab="",main.title="",peak.rate.prefix="")
{
  jet.colors = colorRampPalette(c("#00007F", "blue","#007FFF",  "cyan", "#7FFF7F", "yellow", "#FF7F00","red"))
  par(oma=outma,mar=margin)

  if(class(m)=="matrix"){
    maxRate=max(m,na.rm=T)
    if(maxRate==-1.0)
      image(m,zlim=c(-2,-1), col="white",xlab='',ylab='',axes=FALSE)
    else if (maxRate!=0)
      image(m,zlim=c(0,maxRate), col=jet.colors(200),xlab='',ylab='',axes=FALSE)
    else
      image(m,zlim=c(0,1), col=jet.colors(200),xlab='',ylab='',axes=FALSE)
    mtext(paste(peak.rate.prefix,name,round(max(m,na.rm=T),digits=2),"Hz"),line=-0.15,cex=cex.name,side=3)
  }
  if(class(m)=="data.frame"){
    df<-m
    xlen <- length(unique(df$x))
    ylen <- length(unique(df$y))
    zzz <- matrix(df$rate,nrow=xlen,ncol=ylen)
    image(unique(df$x),unique(df$y),zzz,zlim=c(0,max(df$rate,na.rm=T)), col=jet.colors(200),xlab='',ylab='',axes=FALSE)
    mtext(paste(peak.rate.prefix,name,round(max(df$rate,na.rm=T),digits=2),"Hz"),side=3,at=median(unique(df$x)),line=-0.1,cex=cex.name)
  }
  if(main.title!="")
  {
    mtext(main.title,side=3,line=0.35,cex=cex.title)
  }
}


#' Plot several firing rate maps on the same page
#' 
#' Plot 2-dimensional representation of firing rate for several maps
#' This is not currently being developed.
#' 
#' @param maps A 3d array containing maps (x,y,clu)
#' @param names A character vector containing the name of each map in the array
#' @param ncol Number of columns of plots per page
#' @param nrow Number of rows of plots per page
firingRateMapsPlot<-function(maps,names,ncol=5,nrow=6){
  plot.per.page=ncol*nrow
  m<-matrix(c(rep(seq(0,1-(1/ncol),1/ncol),nrow),
              rep(seq(1/ncol,1,1/ncol),nrow),
              rep(seq(1-(1/nrow),0,0-1/nrow),each=ncol),
              rep(seq(1,1/nrow,0-1/nrow),each=ncol)),ncol=4)
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
}




#' Plot a single spatial autocorrelation map
#' 
#' 
#' @param m A matrix containing the data of the spatial autocorrelation.
#' @param name Character vector giving the name of the map
#' @param outma Outer margins of the figure
#' @param margin Inner margins of the figure
#' @param main.title A title for the figure
#' @param peak.rate.prefix Additional information to display before the peak value.
firingRateMapAutoPlot <- function(m,name="",
                              outma=c(2.0,2.0,2.0,2.0),margin=c(1,1,1,1),
                              main.title="",peak.rate.prefix="")
{
  jet.colors = colorRampPalette(c("#00007F", "blue","#007FFF",  "cyan", "#7FFF7F", "yellow", "#FF7F00","red"))
  par(oma=outma,mar=margin)
  image(m,zlim=c(-1,1), col=jet.colors(200),xlab='',ylab='',axes=FALSE)
  mtext(paste(peak.rate.prefix,name,round(max(m,na.rm=T),digits=2),"r value"),line=-0.1,cex=0.6,side=3)
  if(main.title!="")
  {
    mtext(main.title,side=3,line=0.3,cex=0.5)
  }
}


#' Plot a several spatial autocorrelation maps on the same page
#' 
#' This is not currently being developed.
#' 
#' @param maps A 3d array containing maps (x,y,clu)
#' @param names A character vector containing the name of each map in the array
#' @param ncol Number of columns of plots per page
#' @param nrow Number of rows of plots per page
firingRateMapAutosPlot<-function(maps,names,ncol=5,nrow=6){
  plot.per.page=ncol*nrow
  m<-matrix(c(rep(seq(0,1-(1/ncol),1/ncol),nrow),
              rep(seq(1/ncol,1,1/ncol),nrow),
              rep(seq(1-(1/nrow),0,0-1/nrow),each=ncol),
              rep(seq(1,1/nrow,0-1/nrow),each=ncol)),ncol=4)
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
}



#' Plot a single spatial crosscorrelation map
#' 
#' 
#' @param m A matrix containing the data of the spatial crosscorrelation.
#' @param name Character vector giving the name of the map
#' @param outma Outer margins of the figure
#' @param margin Inner margins of the figure
#' @param main.title A title for the figure
#' @param peak.rate.prefix Additional information to display before the peak value.
firingRateMapCrossPlot <- function(m,name="",
                                  outma=c(2.0,2.0,2.0,2.0),margin=c(1,1,1,1),
                                  main.title="",peak.rate.prefix="")
{
  jet.colors = colorRampPalette(c("#00007F", "blue","#007FFF",  "cyan", "#7FFF7F", "yellow", "#FF7F00","red"))
  par(oma=outma,mar=margin)
  image(m,zlim=c(-1,1), col=jet.colors(200),xlab='',ylab='',axes=FALSE)
  graphics::segments(0.5,0.0,0.5,1)
  graphics::segments(0.0,0.5,1.0,0.5)
  mtext(paste(peak.rate.prefix,name,round(max(m,na.rm=T),digits=2),"r value"),line=-0.1,cex=0.6,side=3)
  if(main.title!="")
  {
    mtext(main.title,side=3,line=0.3,cex=0.5)
  }
}

















#' Polar plot of firing rate as a function of head direction
#' 
#' Polar representation of firing rate using the polar.plot function of the plotrix package
#' 
#' @param histo A data.frame with columns deg and rate, or a numeric vector containing rate values
#' @param outma Outer margins of the figure
#' @param margin Inner margins of the figure
#' @param axis.x.mgp mgp for x axis
#' @param axis.y.mgp mgp for y axis
#' @param cex.x.axis cex for x axis
#' @param cex.y.axis cex for y axis
#' @param cex.lab cex for labels
#' @param xlab Name to display under the x axis
#' @param ylab Name to display at the left of the y axis
#' @param show.xlab Wether of not the xlabel is shown
#' @param main.title A title for the figure
#' @param peak.rate.prefix Additional information to display before the peak firing rate.
#' @param clockwise plot degrees clockwise or anticlockwise
headDirectionPolarPlot <- function(histo,outma=c(0.5,0.5,0.5,0.5),margin=c(0.5,0.3,0.5,0.3),axis.x.mgp=c(1,0.3,0),
                                   axis.y.mgp=c(2.2,0.6,0),cex.x.axis=0.5,cex.y.axis=0.5,cex.lab=0.5,
                                   xlab="",ylab="",show.xlab=TRUE,main.title="",peak.rate.prefix="",
                                   clockwise=TRUE)
{
  par(oma=outma,
      mar=margin,
      cex.lab=cex.lab,cex.axis=cex.x.axis,cex.lab=cex.lab)
  
  if(class(histo)=="numeric"){
    int<-360/length(histo)
    deg<-seq(0,360-int,int)
    radlim=max(histo)   
    oldpar <- plotrix::polar.plot(histo,
                                  polar.pos=deg,
                                  labels=seq(0,270,90),label.pos=c(0,90,180,270),start=0,
                                  #                       labels="",
                                  clockwise=clockwise,
                                  rp.type="p",
                                  show.grid=T,show.radial.grid=T,
                                  radial.lim=c(0,radlim),show.grid.labels=0,
                                  xlab="",ylab="",line.col=4,mar=margin)
  }
  if(class(histo)=="data.frame"){
    if(!"rate" %in% names(histo))
      stop("headDirectionPolarPlot needs a data.frame with rate column")
    if(!"deg" %in% names(histo))
      stop("headDirectionPolarPlot needs a data.frame with deg column")
    radlim=max(histo$rate)   
    oldpar <- plotrix::polar.plot(histo$rate,
                                  polar.pos=histo$deg,
                                  labels=seq(0,270,90),label.pos=c(0,90,180,270),start=0,
                                  #                       labels="",
                                  clockwise=clockwise,
                                  rp.type="p",
                                  show.grid=T,show.radial.grid=T,
                                  radial.lim=c(0,radlim),show.grid.labels=0,
                                  xlab="",ylab="",line.col=4,mar=margin)
  }
  mtext(paste(peak.rate.prefix,round(radlim,digits=2),"Hz"),side=3,at=0,line=0.1,cex=0.5)  
  if(main.title!="")
    mtext(main.title,side=3,at=0,line=0.3,cex=0.5)
  par(oldpar)
}


#' Plot several polar plot with firing rate as a function of head direction
#' 
#' 
#' @param histo A matrix [deg,clu], or numeric if there is only one cell, containing the data for the polar plots.
#' @param names A character vector containing the name of each plot
#' @param ncol Number of columns of plots per page
#' @param nrow Number of rows of plots per page
headDirectionPolarPlots<-function(histo,names,ncol=5,nrow=6){
  plot.per.page=ncol*nrow
  m<-matrix(c(rep(seq(0,1-(1/ncol),1/ncol),nrow),
              rep(seq(1/ncol,1,1/ncol),nrow),
              rep(seq(1-(1/nrow),0,0-1/nrow),each=ncol),
              rep(seq(1,1/nrow,0-1/nrow),each=ncol)),ncol=4)
  index=1
  if(class(histo)=="numeric")
    nCells<-1
  if(class(histo)=="matrix")
    nCells<-dim(histo)[2]
  
  for (i in 1:nCells){
    if(index==1)
    {
      split.screen(m)  
    }
    screen(index)
    ## insert your plot function here
    if(nCells==1)
      headDirectionPolarPlot(histo,peak.rate.prefix=names[i])
    else
      headDirectionPolarPlot(histo[,i],peak.rate.prefix=names[i])
    if(index==plot.per.page)
    {
      close.screen( all.screens = TRUE )
      index=0
    }
    index=index+1
  }
  close.screen(all.screens = TRUE)
}

#' Plot a linear rate histogram from a SpatialProperties1D object
#' 
#' In development
#' 
#' @param sp1 SpatialProperties1d object
#' @param n Index of the histogram in the sp1 object that you want to plot
#' @param outma Outer margins of the figure
#' @param margin Inner margins of the figure
#' @param axis.x.pos position of x axis
#' @param axis.y.pos position of y axis
#' @param axis.y.las las of y axis
#' @param mgp.x mgp for x axis
#' @param mgp.y mgp for y axis
#' @param xlab Name to display under the x axis
#' @param ylab Name to display at the left of the y axis
#' @param xaxis.at where to put the tics
#' @param yaxis.at where to put the tics
#' @param ... passed to the plot function
linearRatePlot<-function(sp1,n=1,
                         outma=c(1,1,1,0),margin=c(1.5,2,0.8,0.3),
                         axis.x.pos=0,axis.y.pos=0,axis.y.las=2,
                           mgp.x=c(0.5,0.05,0.0),mgp.y=c(.8,0.3,0.2),xlab="Position (cm)",ylab="Rate (Hz)",
                           xaxis.at=seq(0,80,20),yaxis.at=seq(0,50,10),...)
{
  rate=sp1@rateHisto[,n]
  rate[which(rate==-1.0)]<-NA
  position = seq(sp1@cmPerBin,sp1@cmPerBin*sp1@nBinRateHisto,sp1@cmPerBin)
  plotxlim=c(0,max(position))
  plotylim=c(0,max(rate,na.rm=T))
  par(mar=margin, oma=outma,cex.lab=0.6,cex.axis=0.6)
  plot (x=plotxlim, y=plotylim,type='n', axes=FALSE, pch=20,lwd=1,xlab="",ylab="",...)
  axis(side = 1, pos=axis.x.pos, at=xaxis.at, tck=-0.05,cex.axis=0.60,mgp=mgp.x)
  par(mgp=mgp.y)
  axis(side = 2, las=axis.y.las, pos=axis.y.pos,tck=-0.05,cex.axis=0.60,mgp=mgp.y)
  lines (position,rate,type='l',pch=20,xlab='',ylab='',lwd=0.75,col="black")
  title(xlab=xlab,mgp=mgp.x)
  title(ylab=ylab,mgp=mgp.y)
}

#' Plot several linear rate histograms from a SpatialProperties1D object
#' 
#' In development
#' 
#' @param sp1 SpatialProperties1d object
#' @param ncol Number of columns of plots per page
#' @param nrow Number of rows of plots per page
linearRatePlots<-function(sp1,ncol=5,nrow=6){
  plot.per.page=ncol*nrow
  m<-matrix(c(rep(seq(0,1-(1/ncol),1/ncol),nrow),
              rep(seq(1/ncol,1,1/ncol),nrow),
              rep(seq(1-(1/nrow),0,0-1/nrow),each=ncol),
              rep(seq(1,1/nrow,0-1/nrow),each=ncol)),ncol=4)
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
}
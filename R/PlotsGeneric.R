#' Plot a boxplot with 2 factors
#' 
#' Just a wrapper around boxplot function.
#' To set which level is plot first, consider using the function factor and set the argument levels.
#' 
#' @param data data.frame with the data
#' @param dv Dependent variable. The one you want on the y-axis
#' @param iv1 First independent variable
#' @param iv2 Second independent variable
#' @param colors Colors of the bars of the boxplot
#' @param ylim Limit of y axis
#' @param ylab Label of y axis
#' @param outma outma for setting par
#' @param margin margin for setting par
#' @param mgp mgp for setting par
#' @param legend.text Text for the legend
#' @param legend.color List of color for the legend
#' @param legend.xy x and y position of the legend
#' @param ... passed to graphics::boxplot function
boxplotTwoFactors <- function(data,
                              dv="vd",
                              iv1="vi1",
                              iv2="vi2",
                              colors=c("red","gray","red","gray"),
                              ylim=c(0,1),
                              ylab="",
                              outma=c(0.5,1,0.5,1),
                              margin=c(1,2,1,1),
                              mgp=c(1,0.2,0),
                              legend.text="",
                              legend.color="",
                              legend.xy="",
                              ...)
{
  par(oma=outma,mar=margin,mgp=mgp)
  graphics::boxplot(data[,dv]~data[,iv1]*data[,iv2],las=1,
          frame.plot=FALSE,axes=FALSE,outline=FALSE,
          col=colors,
          at=c(1,2,4,5),ylim=ylim,...)
  # axis(side = 1, at=0:2,pos=0, tck=-0.05,cex.axis=0.6,labels=c("","",""))
  par(mgp=mgp)
  graphics::axis(side = 2,las=1, pos=0,tck=-0.05,cex.axis=0.6,xpd=TRUE)
  graphics::title(ylab=ylab,mgp=mgp,cex=0.6)
  if(legend.text[1]!=""){
    par(xpd=T) # remove clipping so part of the legend can be outside the plot
    graphics::legend(legend.xy[1], legend.xy[2], legend=legend.text,col=legend.color, lty=1, lwd = 2, cex=0.6,bty = "n")
  }
}

#' Plot points in 2d
#' 
#' Just a wrapper around plot and points functions.
#' 
#' @param data data.frame with the data
#' @param v1 Name of the variable to plot on the x axis
#' @param v2 Name of the variable to plot on the y axis
#' @param axis.y.pos Position of the y axis
#' @param axis.x.pos Position of the x axis
#' @param axis.y.las Orientation of the letters for the label of the y axis
#' @param main.title Main title
#' @param mgp.x mgp for x axis
#' @param mgp.y mgp for y axis
#' @param xlab Label for the x axis
#' @param ylab Label for the y axis
#' @param plotxlim Limit of the x axis
#' @param plotylim Limit of the y axis
#' @param outma outma for setting par
#' @param margin margin for setting par
#' @param xaxis.at Where the tics are shown on the x axis
#' @param yaxis.at Where the tics are shown on the y axis
#' @param cex.point cex for the size of the points
#' @param add.text Text to add to the plot
#' @param add.text.pos Position of the text to add, format (x,y)
#' @param ... Passed to the graphics::plot function
plotPoints <- function(data,v1="v1",v2="v2",axis.y.pos=-.2,axis.x.pos=-.2,axis.y.las=2,
                        main.title="",mgp.x=c(0.5,0.1,0.1),mgp.y=c(0.9,0.2,0.1),xlab="",ylab="",
                        plotxlim=c(-.2,0.6),plotylim=c(-0.2,0.6),outma=c(0.5,0.5,0.5,0.5),margin=c(1.5,1.5,1,0.3),
                        xaxis.at=seq(-0.2,0.7,.2),yaxis.at=seq(-0.2,0.7,.2),cex.point=0.1,
                        add.text="",
                        add.text.pos=c(0,0.5),...)
{
  par(mar=margin, oma=outma,cex.lab=0.6,cex.axis=0.6)
  graphics::plot(x=plotxlim,y=plotylim,type='n', axes=FALSE, pch=20,lwd=1,xlab="",ylab="",...)
  points(data[,v1],data[,v2],pch=20,cex=cex.point)
  par(mgp=mgp.x)
  graphics::axis(side = 1, at=xaxis.at, pos=axis.x.pos,tck=-0.05,cex.axis=0.6)
  par(mgp=mgp.y)
  graphics::axis(side = 2, at=yaxis.at, las=axis.y.las, pos=axis.y.pos,tck=-0.05,cex.axis=0.6)
  graphics::title(xlab=xlab,mgp=mgp.x)
  graphics::title(ylab=ylab,mgp=mgp.y)
  if(main.title!=""){
    graphics::title(main=main.title,cex.main=0.4)
  }
  if(add.text!=""){
    graphics::text(labels=add.text,x=add.text.pos[1],y=add.text.pos[2],cex=0.6)
  }
}




#' Plot a single spike-time autocorrelation
#' 
#' @param y Numeric vectors with the y values
#' @param x Numeric vectors with the x values
#' @param name Character vectors containing the name of the graph
#' @param axis.y.pos Position of the y axis
#' @param axis.x.pos Position of the x axis
#' @param axis.y.las Orientation of the letters for the label of the y axis
#' @param main.title Main title
#' @param mgp.x mgp for x axis
#' @param mgp.y mgp for y axis
#' @param xlab Label for the x axis
#' @param ylab Label for the y axis
#' @param plotxlim Limit of the x axis
#' @param plotylim Limit of the y axis
#' @param outma outma for setting par
#' @param margin margin for setting par
#' @param xaxis.at Where the tics are shown on the x axis
#' @param yaxis.at Where the tics are shown on the y axis
#' @param add.text Text to add to the plot
#' @param add.text.pos Position of the text to add, format (x,y)
#' @param ... Passed to the graphics::plot function
spikeTimeAutocorrelationPlot <- function(y,x,name="",
                                         axis.y.pos=NA,axis.x.pos=0,
                                         axis.y.las=2,
                                         main.title="",mgp.x=c(0.5,0.1,0.1), mgp.y=c(0.9,0.2,0.1),
                                         xlab="Time (ms)", ylab="Spikes",
                                         plotxlim=NA,plotylim=NA,
                                         outma=c(0.5,0.5,0.5,0.5),margin=c(1.5,1.5,1,0.3),
                                         xaxis.at=NA,yaxis.at=NA,
                                         add.text="",add.text.pos=c(0,0.5),...)
{
  
   par(mar=margin, oma=outma,cex.lab=0.6,cex.axis=0.6)
  if(is.na(plotxlim))
      plotxlim=c(min(x),max(x))
  if(is.na(plotylim))
    plotylim=c(0,max(y))
  if(is.na(axis.y.pos))
    axis.y.pos<-min(x)
  if(length(x)!=length(y))
    stop("length(x) != length(y)")
  graphics::plot(x=plotxlim,y=plotylim,type='n', axes=FALSE, pch=20,lwd=1,xlab="",ylab="") #,...)
  lines(x,y)
  par(mgp=mgp.x)
  if(is.na(xaxis.at)){
    graphics::axis(side = 1, pos=axis.x.pos, tck=-0.05,cex.axis=0.6)
  }else{
    graphics::axis(side = 1, pos=axis.x.pos, at=xaxis.at, tck=-0.05,cex.axis=0.6)
  }
  par(mgp=mgp.y)
  if(is.na(yaxis.at)){
    graphics::axis(side = 2, las=axis.y.las, pos=axis.y.pos,tck=-0.05,cex.axis=0.6)
  } else{
    graphics::axis(side = 2, at=yaxis.at, las=axis.y.las, pos=axis.y.pos,tck=-0.05,cex.axis=0.6)
  }
  
  graphics::title(xlab=xlab,mgp=mgp.x)
  graphics::title(ylab=ylab,mgp=mgp.y)
  if(main.title!=""){
    graphics::title(main=main.title,cex.main=0.4)
  }
  if(add.text!=""){
    graphics::text(labels=add.text,x=add.text.pos[1],y=add.text.pos[2],cex=0.6)
  }
}

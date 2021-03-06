#' A S4 class to represent the spike waveforms of neurons
#'
#' This class is used to get the mean waveform of a neuron. 
#' It can also extract some features from the spike waveform.
#'
#'
#' @slot session Name of the recording session
#' @slot path Directory where the recording session is located
#' @slot samplingRate Sampling rate of the electrophysiological data
#' @slot cellList Cell list
#' @slot wf Matrix holding spike waveforms of the neurons
#' @slot wfCluId Character vector indicating which waveforms belong to which neuron
#' @slot wfMsPerBin Ms per bin in spike-time autocorrelation
#' @slot wfTimePoints Time points for data points in the spike-time autocorrelation

SpikeWaveform <- setClass(
  "SpikeWaveform", ## name of the class
  slots=c(session="character",
          path="character",
          samplingRate="numeric",
          cellList="numeric",
          wf="matrix",
          wfCluId="character",
          wfMsPerBin="numeric",
          wfTimePoints="numeric"
          ),
  prototype = list(session=""))


#' Calculate the mean waveform of neurons in a recording session.
#'
#' @param sw SpikeWaveform object
#' @param rs RecSession object
#' @param st SpikeTrain object
#' @param cg CellGroup object
#' @param df DatFiles object
#' @param windowSizeMs The window size for the mean waveform
#' @param divisorToMicroVolt Factor by which to divide the raw data to obtain values in micro volts.
#' @param numberSpikes Number of spikes to use to calculate the mean
#' @param minInterSpikeIntervalMs Minimum inter spike intervals between the spikes used to calculate the mean
#' @param minPassHz Minimal frequency of the bandpass filter, only used if filter==TRUE
#' @param maxPassHz Maximal frequency of the bandpass filter, only used if filter==TRUE
#' @param filter Logical indicating whether to apply a bandpass filter to raw traces
#' per neuron so that all neurons have the same number of channels.

#' @return SpikeWaveform object with the mean waveform in wf matrix
#' 
#' @docType methods
#' @rdname meanWaveform-methods
setGeneric(name="meanWaveform",
           def=function(sw,rs,st,cg,df,windowSizeMs=3,divisorToMicroVolt=7.0,
                        numberSpikes=500,minInterSpikeIntervalMs=2,minPassHz=500,maxPassHz=10000,filter=FALSE)
           {standardGeneric("meanWaveform")})

#' @rdname meanWaveform-methods
#' @aliases meanWaveform,ANY,ANY-method
setMethod(f="meanWaveform",
          signature = "SpikeWaveform",
          definition=function(sw,rs,st,cg,df,windowSizeMs=3,divisorToMicroVolt=7.0,
                              numberSpikes=500,minInterSpikeIntervalMs=2,minPassHz=500,maxPassHz=10000,filter=FALSE)
          {
            if(rs@session=="")
            {
              stop("RecSession object not set in meanWaveform()")
            }
            if(st@session=="")
            {
              stop("SpikeTrain object not set in meanWaveform()")
            }
            if(windowSizeMs<1000/rs@samplingRate)
            {
              stop(paste("windowSizeMs of",windowSizeMs,
                         "is smaller than one data point given the sampling rate of", 
                         rs@samplingRate,"in meanWaveform()"))
            }
            if(numberSpikes<1)
            {
              stop(paste("numberSpikes",numberSpikes,"is smaller than 1 in meanWavefrom()"))
            }
            if(st@nCells<1)
            {
              stop(paste("st@nCells is smaller than 1 in meanWaveform()"))
            }
            if(length(df@fileNames)==0)
            {
              stop(paste("df@fileNames has a length of 0 in meanWaveform"))
            }
            
            # set some variables of the SpikeWaveform object
            sw@samplingRate=rs@samplingRate
            sw@path=rs@path
            sw@session=rs@session
            sw@cellList=st@cellList
            sw@wfCluId<-character()
            
            
            for(clu in st@cellList)
            {
              #print(clu)
              ## get the spike times
              spikeTimes<-st@res[which(st@clu==clu)]
              
              ## only keep spikes with a valid inter spike interval
              spikeTimes<-spikeTimes[(c(FALSE,diff(spikeTimes)>minInterSpikeIntervalMs*sw@samplingRate/1000))]
              
              ## only keep the number of spikes needed
              spikeTimes<-head(spikeTimes,n=numberSpikes)
              
              ## get the channels on which the cluter was recorded
              tetrode<-cg@tetrode[which(cg@clu==clu)]
              channels<-rs@channelsTetrode[tetrode,]
              channels<-channels[which(!is.na(channels))]
              ## get the raw data for the 4 channels from 0 to last spikes + have time window
              startIndex<-0
              endIndex=tail(spikeTimes,n=1)+(windowSizeMs*rs@samplingRate/1000)
            
              ## load the raw trace from dat files
              traces<-datFilesGetChannels(df,channels,firstSample=startIndex,lastSample=endIndex)
              
              ## filter here so that we can extract the waveform from it later
              if(filter==TRUE)
                for(chan in 1:length(channels))
                  traces[,chan]<-bandPassFilter(as.numeric(traces[,chan]),rs@samplingRate,minPassHz,maxPassHz)
              wf<-spikeWaveformFromTraces(traces,spikeTimes,as.integer(windowSizeMs*rs@samplingRate/1000))
              mwf<-apply(wf,c(2,3),mean)
              if(length(sw@wfCluId)==0)
              {
                sw@wf<-mwf
                sw@wfCluId<-rep(cg@id[which(cg@clu==clu)],length(channels))
              }else{
                sw@wf<-cbind(sw@wf,mwf)
                sw@wfCluId<-c(sw@wfCluId,rep(cg@id[which(cg@clu==clu)],length(channels)))
              }
            }
            sw@wfMsPerBin<-1000/rs@samplingRate
            sw@wfTimePoints<- seq(from = -sw@wfMsPerBin*dim(sw@wf)[1]/2+sw@wfMsPerBin/2, 
                                  to =sw@wfMsPerBin*dim(sw@wf)[1]/2-sw@wfMsPerBin/2, 
                                  by =sw@wfMsPerBin)
            rownames(sw@wf)<-sw@wfTimePoints
           
            if(length(sw@wfCluId)!=dim(sw@wf)[2])
              stop("problem with the dimensions of sw@wf and length of sw@wfCluId in meanWavefrom()")
            
            if(length(sw@wfTimePoints)!=dim(sw@wf)[1])
              stop("problem with the dimensions of sw@wf and length of sw@wfTimePoints in meanWavefrom()")
            return(sw)
          }
)

#' Calculate waveform characteristics from the mean waveform of each cluster.
#'
#' You first need to run meanWavefrom() before calling this function.
#' The wire with the largest range is used.
#' The characteristics that are calculated are amplitude, duration,
#' duration before trough, duration after trough,
#' peak from baseline before trough,
#' peak from baseline after trough,
#' peak amplitude assymetry
#'
#' @param sw SpikeWaveform object
#' @return data.frame with the different waveform features for each neuron
#' 
#' @docType methods
#' @rdname waveformCharacteristics-methods
setGeneric(name="waveformCharacteristics",
           def=function(sw)
           {standardGeneric("waveformCharacteristics")})

#' @rdname waveformCharacteristics-methods
#' @aliases waveformCharacteristics,ANY,ANY-method
setMethod(f="waveformCharacteristics",
          signature = "SpikeWaveform",
          definition=function(sw)
          {
            if(dim(sw@wf)[1]==0)
            {
              stop("you need to call meanWaveform() before calling waveformCharacteristics()")
            }
            
            df<-data.frame()
            for(id in unique(sw@wfCluId))
            {
              ## get the waveforms for this CluId
              m<-sw@wf[,which(sw@wfCluId==id)]
              ## select waveform with most negative values
              wf<-m[,which.min(apply(m,2,min))]
              df<-rbind(df,spikeGeoFeatureOneSpike(wf,sw@wfTimePoints))
            }
            df$cluId<-unique(sw@wfCluId)
            colnames(df)<-c("baseline","amplitude","spikeDuration", "firstHalfSpikeDuration", "secondHalfSpikeDuration","peakAmplitudeAssymetry","cluId")
        
            return(df)
          }
)

#' Plot the mean waveform of a cluster using a spikeWaveform object
#'
#' You first need to run meanWavefrom() before calling this function.
#'
#' @param sw SpikeWaveform object
#' @param cluId cluId of the mean waveform. cluId contains the fileBase and cluster number seperated by underscore
#' @param offset offset between the traces in the plot
#' @param name name for the plot
#' @param axis.y.pos position of the y axis
#' @param axis.x.pos position of the x axis
#' @param axis.y.las orientation of numbers on y axis, 1 or 2
#' @param mgp.x mgp of x axis
#' @param mgp.y mgp of y axis
#' @param xlab label of x axis
#' @param ylab label of y axis
#' @param plotxlim limits of x axis
#' @param plotylim limits of y axis
#' @param outma outma for graph
#' @param margin margin of graph
#' @param xaxis.at x axis tic location
#' @param yaxis.at y axis tics location
#' @param add.text text you can add on the graph
#' @param add.text.pos position of the text you can add on the graph
#' 
#' @docType methods
#' @rdname spikeWaveformPlot-methods
setGeneric(name="spikeWaveformPlot",
           def=function(sw,cluId,offset=200,
                        name="",
                        axis.y.pos=NA,
                        axis.x.pos=NA,
                        axis.y.las=2,
                        mgp.x=c(0.5,0.1,0.1),
                        mgp.y=c(1.1,0.2,0.1),
                        xlab="Time (ms)",
                        ylab="Voltage",
                        plotxlim=NA,
                        plotylim=NA,
                        outma=c(0.5,0.5,0.5,0.5),
                        margin=c(1.5,1.7,1,0.3),
                        xaxis.at=NA,
                        yaxis.at=NA,
                        add.text="",
                        add.text.pos=c(0,0.5))
           {standardGeneric("spikeWaveformPlot")})

#' @rdname spikeWaveformPlot-methods
#' @aliases spikeWaveformPlot,ANY,ANY-method
setMethod(f="spikeWaveformPlot",
          signature = "SpikeWaveform",
          definition=function(sw,cluId,offset=250,
                              name="",
                              axis.y.pos=NA,
                              axis.x.pos=NA,
                              axis.y.las=2,
                              mgp.x=c(0.5,0.1,0.1),
                              mgp.y=c(1.1,0.2,0.1),
                              xlab="Time (ms)",
                              ylab="Voltage",
                              plotxlim=NA,
                              plotylim=NA,
                              outma=c(0.5,0.5,0.5,0.5),
                              margin=c(1.5,1.7,1,0.3),
                              xaxis.at=NA,
                              yaxis.at=NA,
                              add.text="",
                              add.text.pos=c(0,0.5))
          {
            if(dim(sw@wf)[1]==0)
            {
              stop("you need to call meanWaveform() before calling waveformPlot()")
            }
            
            if(cluId%in%sw@wfCluId==FALSE)
            {
              stop(paste("cluId",cluId, "is not in the spikeWaveform object"))
            }
            ## get the waveforms for this CluId
            m<-sw@wf[,which(sw@wfCluId==cluId)]
            nChannels<-dim(m)[2]
            offset<-offset*(0:(nChannels-1))
            
            par(mar=margin, oma=outma,cex.lab=0.6,cex.axis=0.6)
            if(any(is.na(plotxlim)))
              plotxlim=c(min(sw@wfTimePoints),max(sw@wfTimePoints))
            if(any(is.na(plotylim)))
              plotylim=c(min(m),max(m)+max(offset))
            if(any(is.na(axis.y.pos)))
              axis.y.pos<-min(sw@wfTimePoints)
            if(any(is.na(axis.x.pos)))
              axis.x.pos<-min(m)
            graphics::plot(x=plotxlim,y=plotylim,type='n', axes=FALSE, pch=20,lwd=1,xlab="",ylab="")
            
            #keep only data within xlim
            m<-m[which(sw@wfTimePoints>=plotxlim[1]&sw@wfTimePoints<=plotxlim[2]),]
            t<-sw@wfTimePoints[which(sw@wfTimePoints>=plotxlim[1]&sw@wfTimePoints<=plotxlim[2])]
            
            for(i in 1:nChannels)
              lines(t,m[,i]+offset[i],col=i)
            par(mgp=mgp.x)
            if(is.na(xaxis.at)){
              graphics::axis(side = 1, pos=axis.x.pos, tck=-0.05, cex.axis=0.6)
            }else{
              graphics::axis(side = 1, pos=axis.x.pos, at=xaxis.at, tck=-0.05, cex.axis=0.6)
            }
            par(mgp=mgp.y)
            if(is.na(yaxis.at)){
              graphics::axis(side = 2, las=axis.y.las, pos=axis.y.pos,tck=-0.05,cex.axis=0.6)
            } else{
              graphics::axis(side = 2, at=yaxis.at, las=axis.y.las, pos=axis.y.pos,tck=-0.05,cex.axis=0.6)
            }
            graphics::title(xlab=xlab,mgp=mgp.x)
            graphics::title(ylab=ylab,mgp=mgp.y)
            if(name!=""){
              graphics::title(main=name,cex.main=0.5)
            }
            if(add.text!=""){
              graphics::text(labels=add.text,x=add.text.pos[1],y=add.text.pos[2],cex=0.6)
            }
            return()
          }
)




#' Plot the mean waveform of a cluster using a spikeWaveform object
#'
#'
#' @param m matrix with the spike waveform
#' @param wfTimePoints vector with the time points for the waveforms
#' @param offset offset between the traces in the plot
#' @param name name for the plot
#' @param axis.y.pos position of the y axis
#' @param axis.x.pos position of the x axis
#' @param axis.y.las orientation of numbers on y axis, 1 or 2
#' @param mgp.x mgp of x axis
#' @param mgp.y mgp of y axis
#' @param xlab label of x axis
#' @param ylab label of y axis
#' @param plotxlim limits of x axis
#' @param plotylim limits of y axis
#' @param outma outma for graph
#' @param margin margin of graph
#' @param xaxis.at x axis tic location
#' @param yaxis.at y axis tics location
#' @param add.text text you can add on the graph
#' @param add.text.pos position of the text you can add on the graph
waveformPlot<-function(m,
                       wfTimePoints,
                       offset=250,
                       name="",
                       axis.y.pos=NA,
                       axis.x.pos=NA,
                       axis.y.las=2,
                       mgp.x=c(0.5,0.1,0.1),
                       mgp.y=c(1.1,0.2,0.1),
                       xlab="Time (ms)",
                       ylab="Voltage",
                       plotxlim=NA,
                       plotylim=NA,
                       outma=c(0.5,0.5,0.5,0.5),
                       margin=c(1.5,1.7,1,0.3),
                       xaxis.at=NA,
                       yaxis.at=NA,
                       add.text="",
                       add.text.pos=c(0,0.5))
{
  if(dim(m)[1]==0)
  {
    stop("matrix is empty in waveformPlot()")
  }
  if(length(wfTimePoints)==0)
  {
    stop("wfTimePoints has a length of 0")
  }
  if(dim(m)[1]!=length(wfTimePoints))
  {
    stop("wfTimePoints has a length different than the waveforms")
  }

  nChannels<-dim(m)[2]
  offset<-offset*(0:(nChannels-1))

  par(mar=margin, oma=outma,cex.lab=0.6,cex.axis=0.6)
  if(any(is.na(plotxlim)))
    plotxlim=c(min(wfTimePoints),max(wfTimePoints))
  if(any(is.na(plotylim)))
    plotylim=c(min(m),max(m)+max(offset))
  if(any(is.na(axis.y.pos)))
    axis.y.pos<-min(wfTimePoints)
  if(any(is.na(axis.x.pos)))
    axis.x.pos<-min(m)
  graphics::plot(x=plotxlim,y=plotylim,type='n', axes=FALSE, pch=20,lwd=1,xlab="",ylab="")

  #keep only data within xlim
  m<-m[which(wfTimePoints>=plotxlim[1]&wfTimePoints<=plotxlim[2]),]
  t<-wfTimePoints[which(wfTimePoints>=plotxlim[1]&wfTimePoints<=plotxlim[2])]

  for(i in 1:nChannels)
    lines(t,m[,i]+offset[i],col=i)
  par(mgp=mgp.x)
  if(is.na(xaxis.at)){
    graphics::axis(side = 1, pos=axis.x.pos, tck=-0.05, cex.axis=0.6)
  }else{
    graphics::axis(side = 1, pos=axis.x.pos, at=xaxis.at, tck=-0.05, cex.axis=0.6)
  }
  par(mgp=mgp.y)
  if(is.na(yaxis.at)){
    graphics::axis(side = 2, las=axis.y.las, pos=axis.y.pos,tck=-0.05,cex.axis=0.6)
  } else{
    graphics::axis(side = 2, at=yaxis.at, las=axis.y.las, pos=axis.y.pos,tck=-0.05,cex.axis=0.6)
  }
  graphics::title(xlab=xlab,mgp=mgp.x)
  graphics::title(ylab=ylab,mgp=mgp.y)
  if(name!=""){
    graphics::title(main=name,cex.main=0.5)
  }
  if(add.text!=""){
    graphics::text(labels=add.text,x=add.text.pos[1],y=add.text.pos[2],cex=0.6)
  }
  return()
}









#' Get the geometrical features of a single waveform
#' 
#' @param wf Numeric vector with waveform
#' @param timePoints Numeric vector with time points
#' @param baselineProportion Numeric indicating the percentage of wf 
#' @param thresholdAmplitudeProportion
#' that is used to calculate the baseline
#' @return spike geometrical features (baseline,amplitude,spikeDuration, firstHalfSpikeDuration, secondHalfSpikeDuration,peakAmplitudeAssymetry)
spikeGeoFeatureOneSpike<-function(wf,timePoints,baselineProportion=0.15,thresholdAmplitudeProportion=0.5){

  baseline<-mean(wf[1:(length(wf)*baselineProportion)])
  trough=min(wf)
  troughIndex=which.min(wf)
  amplitude=baseline-trough
  thresholdAmplitude=baseline-(amplitude*thresholdAmplitudeProportion)
  startSpikeIndex<-utils::tail(which(wf[1:troughIndex]>thresholdAmplitude),n=1)
  endSpikeIndex<-utils::head(which(wf[troughIndex:length(wf)]>thresholdAmplitude),n=1)+troughIndex-1
  interpolationStart <- (wf[startSpikeIndex]-thresholdAmplitude)/(wf[startSpikeIndex]-wf[startSpikeIndex+1])
  interpolationEnd<-(wf[endSpikeIndex]-thresholdAmplitude)/(wf[endSpikeIndex]-wf[endSpikeIndex-1])
  interpolationStart<-startSpikeIndex+interpolationStart
  interpolationEnd<- endSpikeIndex-interpolationEnd
  spikeDuration=(interpolationEnd-interpolationStart)*(timePoints[2]-timePoints[1])
  firstHalfSpikeDuration=(troughIndex-interpolationStart)*(timePoints[2]-timePoints[1])
  secondHalfSpikeDuration=(interpolationEnd-troughIndex)*(timePoints[2]-timePoints[1])
  maxPreTrough<-max(wf[1:troughIndex])-baseline
  maxPostTrough<-max(wf[troughIndex:length(wf)])-baseline
  peakAmplitudeAssymetry<-(maxPreTrough - maxPostTrough)/(abs(maxPreTrough)+abs(maxPostTrough))
  return(c(baseline,amplitude,spikeDuration, firstHalfSpikeDuration, secondHalfSpikeDuration,peakAmplitudeAssymetry))
}











### show ###
setMethod("show", "SpikeWaveform",
          function(object){
            print(paste("session:",object@session))
            print(paste("path:",object@path))
            print(paste("samplingRate:",object@samplingRate))
            print("cellList:")
            print(object@cellList)
            if(dim(object@wf)[1]!=0){
              print(paste("number of waveforms:",dim(object@wf)[2]))
            }
          })


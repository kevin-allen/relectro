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
#' @slot wfMsPerBin Ms per bin in spike-time autocorrelation
#' @slot wfTimePoints Time points for data points in the spike-time autocorrelation

spikeWaveform <- setClass(
  "SpikeWaveform", ## name of the class
  slots=c(session="character",
          path="character",
          samplingRate="numeric",
          cellList="numeric",
          wf="matrix",
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
#' @minPassHz Minimal frequency of the bandpass filter, only used if filter==TRUE
#' @maxPassHz Maximal frequency of the bandpass filter, only used if filter==TRUE

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
            
            for(clu in st@cellList)
            {
              print(clu)
              
              ## get the spike times
              spikeTimes<-st@res[which(st@clu==clu)]
              
              ## only keep spikes with a valid inter spike interval
              spikeTimes<-spikeTimes[(c(FALSE,diff(spikeTimes)>minInterSpikeIntervalMs*sw@samplingRate/1000))]
              
              ## only keep the number of spikes needed
              spikeTimes<-head(spikeTimes,n=numberSpikes)
              
              ## get the channels on which the cluter was recorded
              tetrode<-cg@tetrode[which(cg@clu==clu)]
              channels<-rs@channelsTetrode[tetrode,]
              
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
              
              
            dim(wf)
            plot(wf[1,,4],type='l')
              
            }
            
            spikeWaveformFromTraces
            spikeGeoFeatures
            
            
            return(sw)
          }
)




### show ###
setMethod("show", "SpikeWaveform",
          function(object){
            print(paste("session:",object@session))
            print(paste("path:",object@path))
            print(paste("samplingRate:",object@samplingRate))
            print("cellList:")
            print(object@cellList)
          })

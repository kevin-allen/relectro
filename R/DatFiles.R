#' An S4 class representing a group of .dat files
#' 
#' This class is used to read data from one or several .dat files.
#' Because the .dat files are binary files, you can't read them with the usual read functions within R.
#' If several dat files are set in an object, they are treated as though they are a single file.
#' @slot fileNames A character vector containing the names of the files
#' @slot path The directory in which the files are located
#' @slot resofs The number of samples in each file.
#' @slot nChannels Number of recorded channels in the dat files
#' @examples
#' df<-new("DatFile")
#' df<-new("DatFiles",fileName="test.dat",path="~/Documents",nChannels=33)
DatFiles <- setClass(
  "DatFiles", ## name of the class
  slots=c(fileNames="character",
          path="character",
          resofs="numeric",
          nChannels="numeric"
          ),
  prototype = list(fileNames="",path="",nChannels=0))



### datFilesSet ###
setGeneric(name="datFilesSet",
           def=function(df,fileNames,path,nChannels)
           {standardGeneric("datFilesSet")}
)
#' Function to set the parameters of a DatFile object
#' @param df A DatFiles object
#' @param fileNames A character vector containing the names of the .dat files
#' @param path A directory where the files are located
#' @param nChannels The number of channels in the dat files
#' @return A DatFiles object with the values set.
setMethod(f="datFilesSet",
          signature="DatFiles",
          definition=function(df,fileNames="",path="",nChannels=0)
          {
            if(length(fileNames)==0)
              stop("fileNames length = 0")
            if(nChannels<=0)
              stop("nChannels should be larger than 0")
            if(path==""&&df@path==""){
              df@path=getwd()
            } 
            
            df@fileNames=fileNames
            df@nChannels=nChannels
            return(df)                        
          })

setGeneric(name="datFilesGetOneChannel",
           def=function(df,channelNo,firstSample,lastSample)
           {standardGeneric("datFilesGetOneChannel")}
)
#' Function to read the data from a group of dat files
#' 
#' 
#' @param df A DatFiles object
#' @param channelNo The channel of interest
#' @param firstSample First sample to retrieve
#' @param lastSample Last sample to retrieve
#' @return A integer vector containing the data read from the dat files.
#' @aliases datFilesGetOneChannel
setMethod(f="datFilesGetOneChannel",
          signature="DatFiles",
          definition=function(df,channelNo,firstSample,lastSample)
          {
            if(length(df@fileNames)==0)
              stop(paste("df@fileNames length = 0"))
            if(df@nChannels<=0)
              stop(paste("df@nChannles should be larger than 0",df@fileNames[1]))
            if(channelNo<0)
              stop(paste("channelNo < 0",df@fileNames[1]))
            if(channelNo>=df@nChannels)
              stop(paste("channelNo >= df@nChannels",df@fileNames[1]))
            if(firstSample<0)
              stop(paste("firstSample<0",df@fileNames[1]))
            if(lastSample<0)
              stop(paste("lastSample<0",df@fileNames[1]))
            if(firstSample>lastSample)
              stop(paste("firstSample > lastSample",df@fileNames[1]))
            if(length(df@fileNames)>100)
              stop(paste("can only works with a maximum of 100 files",df@fileNames[1]))
            if(any(nchar(df@fileNames)>255))
              stop(paste("max file name size is 255",df@fileNames[1]))
            
           
            results<-.Call("group_data_file_si_get_one_channel_cwrap",
                paste(df@path,df@fileNames,sep="/"),
                df@nChannels,
                channelNo,
                firstSample,
                lastSample)
          return(results)
          })


### show ###
setMethod("show", "DatFiles",
          function(object){
            print(paste("path:",object@path))
            print(paste("fileNames:"))
            print(paste(object@fileNames))
            print(paste("nChannels:",object@nChannels))
          })
#' An S4 class representing a group of .dat files
#' 
#' This class is used to read data from one or several .dat files.
#' Because the .dat files are binary files, you can't read them with the usual read functions within R.
#' If several dat files are set in an object, they are treated as though they are a single file.
#' @slot fileNames A character vector containing the names of the files
#' @slot path The directory in which the files are located
#' @slot samples Numeric vector with the number of samples in each file.
#' @slot size Numeric vector with the size of the file in bytes
#' @slot nChannels Number of recorded channels in the dat files
DatFiles <- setClass(
  "DatFiles", ## name of the class
  slots=c(fileNames="character",
          path="character",
          samples="numeric",
          size="numeric",
          nChannels="numeric"
          ),
  prototype = list(fileNames="",path="",nChannels=0))


#' Function to set the parameters of a DatFile object
#' @param df A DatFiles object
#' @param fileNames A character vector containing the names of the .dat files
#' @param path A directory where the files are located
#' @param nChannels The number of channels in the dat files
#' @return A DatFiles object with the values set.
#' @docType methods
#' @rdname datFilesSet-methods
setGeneric(name="datFilesSet",
           def=function(df,fileNames,path,nChannels)
           {standardGeneric("datFilesSet")}
)

#' @rdname datFilesSet-methods
#' @aliases datFilesSet,ANY,ANY-method
#' 
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
            } else {
              df@path<-path
            }
            df@fileNames=fileNames
            df@nChannels=nChannels
            df<-datFilesSamples(df)
            return(df)                        
          })



#' Get the number of samples per file and file size
#' 
#' @param df A DatFiles object
#' @return A DatFiles object with the samples and size set.
#' @docType methods
#' @rdname datFilesSamples-methods
setGeneric(name="datFilesSamples",
           def=function(df)
           {standardGeneric("datFilesSamples")}
)

#' @rdname datFilesSamples-methods
#' @aliases datFilesSamples,ANY,ANY-method
#' 
setMethod(f="datFilesSamples",
          signature="DatFiles",
          definition=function(df)
          {
            if(length(df@fileNames)==0)
              stop("df@fileNames length = 0")
            if(df@nChannels<=0)
              stop("df@nChannels should be larger than 0")
            ## get the number of channels in each file
            df@size<-(file.info(paste(df@path,df@fileNames,sep="/"))$size)
            ## check that the file size can be divided by nChannels
            if(any(df@size%%(df@nChannels*2)!=0)){
              stop(paste("size of a .dat file can't be divided by",df@nChannels*2))
            }
            df@samples<-df@size/(df@nChannels*2)
            return(df)                        
          })

#' Function to read the data from a group of dat files
#' 
#' 
#' @param df A DatFiles object
#' @param channelNo The channel of interest
#' @param firstSample First sample to retrieve
#' @param lastSample Last sample to retrieve
#' @return A integer vector containing the data read from the dat files.
#' @docType methods
#' @rdname datFilesGetOneChannel-methods
setGeneric(name="datFilesGetOneChannel",
           def=function(df,channelNo,firstSample,lastSample)
           {standardGeneric("datFilesGetOneChannel")}
)

#' @rdname datFilesGetOneChannel-methods
#' @aliases datFilesGetOneChannel,ANY,ANY-method
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
            
            df<-datFilesSamples(df)
            if(firstSample>sum(df@samples))
              stop(paste("firstSample is larger than total number of samples:",sum(df@samples),df@fileNames[1]))
            if(lastSample>sum(df@samples))
              stop(paste("lastSample is larger than total number of samples",sum(df@samples),df@fileNames[1]))
            
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
            print(paste("size in bytes:"))
            print(paste(object@size))
            print(paste("samples per file:"))
            print(paste(object@samples))
          })
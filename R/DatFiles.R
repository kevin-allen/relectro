#' A class representing one or more .dat files
#' 
#' This class is used to read data from one or several .dat files.
#' Because the .dat files are binary files, you can't read them with the usual read functions within R.
#' If several dat files are set in an object, they are treated as though they are a single file.
#' 
#' Dat files have a very simple organisation. Each data point has 2 bytes. 
#' There is no header. The data from sample 0 appear first, then sample 1, etc.
#' If c is for channel and s is for sample, a file with 3 channel would look like this
#' s0c0, s0c1, s0c2, s1c0, s1c1, s1c2, s2c0, s2c1, s2c2, etc.
#' 
#' 
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
#' @param path A directory where the files are located. If empty the working directory is used.
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
            if(path==""){
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
            
            ## check that the files exist
            if(any(!file.exists(paste(df@path,df@fileNames,sep="/")))){
              stop(paste("A dat file is missing",
                         paste(df@path,df@fileNames,sep="/")[which(!file.exists(paste(df@path,df@fileNames,sep="/")))],"\n"))
            
            }
            
            ## get the number of channels in each file
            df@size<-(file.info(paste(df@path,df@fileNames,sep="/"))$size)
            
            ## check that the file size can be divided by nChannels
            if(any(df@size%%(df@nChannels*2)!=0)){
              index<-which(df@size%%(df@nChannels*2)!=0)
              stop(paste(paste("size of can't be divided by",df@nChannels*2), 
                         paste(paste(df@path,df@fileNames[index],sep="/"),"size:",df@size[index],"remainder:",df@size[index]%%df@nChannels*2),"\n"))
            }
            df@samples<-df@size/(df@nChannels*2)
            return(df)                        
          })

#' Function to read one channel from a group of dat files
#' 
#' This function returns the data from firstSample to lastSample. 
#' If firstSample and lastSample have a length greater than 1, 
#' the data chunks will be pasted one after the other in the data returned by the function.
#' 
#' @param df A DatFiles object
#' @param channelNo The channel to retrieve the data from. The first channelNo in the file is 0.  
#' @param firstSample Numeric vector containing the first sample to retrieve. 
#' You can give several if you want to get a series of data chunks.
#' First sample in the file is 0
#' @param lastSample Numeric vector containing the last sample to retrieve. You can give several if you want to get a series of data chunks.
#' Index stats at 0
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
            
            
            if(any(file.exists(paste(df@path,df@fileNames,sep="/"))==FALSE)) {
              stop(paste("This file is missing:",df@fileNames[! file.exists(df@fileNames)]," "))
            }
            
            if(length(df@fileNames)==0)
              stop(paste("df@fileNames length = 0"))
            if(df@nChannels<=0)
              stop(paste("df@nChannles should be larger than 0",df@fileNames[1]))
            if(channelNo<0)
              stop(paste("channelNo < 0",df@fileNames[1]))
            if(channelNo>=df@nChannels)
              stop(paste("channelNo >= df@nChannels",df@fileNames[1]))
            if(length(df@fileNames)>100) # from c code
              stop(paste("can only works with a maximum of 100 files",df@fileNames[1]))
            if(any(nchar(df@fileNames)>255)) # from c code
              stop(paste("max file name size is 255",df@fileNames[1]))
            
            df<-datFilesSamples(df)
            
            if(any(firstSample<0))
              stop(paste("firstSample<0",df@fileNames[1]))
            if(any(lastSample<0))
              stop(paste("lastSample<0",df@fileNames[1]))
            if(any(firstSample>lastSample))
              stop(paste("firstSample > lastSample",df@fileNames[1]))
            if(any(firstSample>(sum(df@samples)-1)))
              stop(paste("firstSample(",firstSample,") is larger or equal than total number of samples:",sum(df@samples),df@fileNames[1]))
            if(any(lastSample> (sum(df@samples)-1)))
              stop(paste("lastSample(",lastSample,") is larger or equal than total number of samples",sum(df@samples),df@fileNames[1]))
            
            if(length(firstSample)!=length(lastSample))
              stop(paste("length(firstSample)!=length(lastSample",df@fileNames[1]))
           
            number.samples=sum(lastSample-firstSample+1)
            results<-numeric(length=number.samples)
            startIndex<-1
            for(i in 1:length(firstSample)){
              s<-lastSample[i]-firstSample[i]+1
              results[startIndex:(startIndex+s-1)]<-.Call("group_data_file_si_get_one_channel_cwrap", 
                                                          paste(df@path,df@fileNames,sep="/"),
                                                          df@nChannels, 
                                                          channelNo,
                                                          firstSample[i], 
                                                          lastSample[i])
              startIndex<-startIndex+s
            }
          return(results)
          })

#' Function to read from several channels from a group of dat files
#' 
#' This function returns the data from firstSample to lastSample for a group of channels. 
#' If firstSample and lastSample have a length greater than 1, 
#' the data chunks will be pasted one after the other in the data returned by the function.
#' 
#' @param df A DatFiles object
#' @param channels Numeric indicating which channels are needed. The first channelNo in the file is 0.  
#' @param firstSample Numeric vector containing the first sample to retrieve. 
#' You can give several if you want to get a series of data chunks.
#' First sample in the file is 0
#' @param lastSample Numeric vector containing the last sample to retrieve. You can give several if you want to get a series of data chunks.
#' Index stats at 0
#' @return A matrix containing the data, one channel per column.
#' @docType methods
#' @rdname datFilesGetChannels-methods
setGeneric(name="datFilesGetChannels",
           def=function(df,channels,firstSample,lastSample)
           {standardGeneric("datFilesGetChannels")}
)

#' @rdname datFilesGetChannels-methods
#' @aliases datFilesGetChannels,ANY,ANY-method
setMethod(f="datFilesGetChannels",
          signature="DatFiles",
          definition=function(df,channels,firstSample,lastSample)
          {
            
            if(any(file.exists(paste(df@path,df@fileNames,sep="/"))==FALSE)) {
              stop(paste("This file is missing:",df@fileNames[! file.exists(df@fileNames)]," "))
            }
            
            if(length(df@fileNames)==0)
              stop(paste("df@fileNames length = 0"))
            if(df@nChannels<=0)
              stop(paste("df@nChannles should be larger than 0",df@fileNames[1]))
            if(length(channels)==0)
              stop(paste("Length of channels = 0",df@fileNames[1]))
            if(any(channels<0))
              stop(paste("at least one channel is < 0",df@fileNames[1]))
            if(any(channels>=df@nChannels))
              stop(paste("At least one channel >= df@nChannels",df@fileNames[1]))
            if(length(df@fileNames)>100) # from c code
              stop(paste("can only works with a maximum of 100 files",df@fileNames[1]))
            if(any(nchar(df@fileNames)>255)) # from c code
              stop(paste("max file name size is 255",df@fileNames[1]))
            
            df<-datFilesSamples(df)
            
            if(any(firstSample<0))
              stop(paste("firstSample<0",df@fileNames[1]))
            if(any(lastSample<0))
              stop(paste("lastSample<0",df@fileNames[1]))
            if(any(firstSample>lastSample))
              stop(paste("firstSample > lastSample",df@fileNames[1]))
            if(any(firstSample>(sum(df@samples)-1)))
              stop(paste("firstSample(",firstSample,") is larger or equal than total number of samples:",sum(df@samples),df@fileNames[1]))
            if(any(lastSample> (sum(df@samples)-1)))
              stop(paste("lastSample(",lastSample,") is larger or equal than total number of samples",sum(df@samples),df@fileNames[1]))
            
            if(length(firstSample)!=length(lastSample))
              stop(paste("length(firstSample)!=length(lastSample",df@fileNames[1]))
            
            number.samples=sum(lastSample-firstSample+1)
            results<-matrix(ncol=length(channels), nrow=number.samples)
            startIndex<-1
            for(i in 1:length(firstSample)){
              s<-lastSample[i]-firstSample[i]+1
              results[startIndex:(startIndex+s-1),]<-
                .Call("group_data_file_si_get_group_channels_cwrap", 
                      paste(df@path,df@fileNames,sep="/"),
                      df@nChannels, 
                      as.integer(channels), 
                      length(channels),
                      firstSample[i], 
                      lastSample[i])
              startIndex<-startIndex+s
            }
            return(results)
          })



#' Function to merge the dat files into a single dat file
#' 
#' 
#' @param df A DatFiles object
#' @param fileName Name of merged file
#' @docType methods
#' @rdname datFilesMerge-methods
setGeneric(name="datFilesMerge",
           def=function(df,fileName)
           {standardGeneric("datFilesMerge")}
)
#' @rdname datFilesMerge-methods
#' @aliases datFilesMerge,ANY,ANY-method
setMethod(f="datFilesMerge",
          signature="DatFiles",
          definition=function(df,fileName)
          {
            if(any(file.exists(paste(df@path,df@fileNames,sep="/"))==FALSE)) {
              stop(paste("This file is missing:",df@fileNames[! file.exists(df@fileNames)]," "))
            }
            if(length(df@fileNames)==0)
              stop(paste("df@fileNames length = 0"))
            myCmd<-paste("cat",  
                         paste(paste(df@path,df@fileNames,sep="/"),collapse=" ") , 
                         ">", 
                         fileName)
            print(paste("Running",myCmd))
            system(myCmd)
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



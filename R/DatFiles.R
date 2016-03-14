############################################
#### definition of DatFiles Class      ###
############################################
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
setMethod(f="datFilesGetOneChannel",
          signature="DatFiles",
          definition=function(df,channelNo,firstSample,lastSample)
          {
            if(length(df@fileNames)==0)
              stop("df@fileNames length = 0")
            if(df@nChannels<=0)
              stop("df@nChannles should be larger than 0")
            if(channelNo<0)
              stop("channelNo < 0")
            if(channelNo>=df@nChannels)
              stop("channelNo >= df@nChannels")
            if(firstSample<0)
              stop("firstSample<0")
            if(lastSample<0)
              stop("lastSample<0")
            if(firstSample>lastSample)
              stop("firstSample > lastSample")
            if(length(df@fileNames)>100)
              stop("can only works with a maximum of 100 files")
            if(any(nchar(df@fileNames)>255))
              stop("max file name size is 255")
            
           
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



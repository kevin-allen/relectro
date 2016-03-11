############################################
#### definition of DatFiles Class      ###
############################################
DatFiles <- setClass(
  "DatFiles", ## name of the class
  slots=c(file.names="character",
          resofs="numeric",
          num.channels="numeric"
          ),
  prototype = list(file.names="",num.channels=0))

### dat.files.set ###
setGeneric(name="dat.files.set",
           def=function(df,file.names,num.channels)
           {standardGeneric("dat.files.set")}
)
setMethod(f="dat.files.set",
          signature="DatFiles",
          definition=function(df,file.names="",num.channels=0)
          {
            if(length(file.names)==0)
              stop("file.names length = 0")
            if(num.channels<=0)
              stop("num.channels should be larger than 0")
            
            df@file.names=file.names
            df@num.channels=num.channels
            return(df)                        
          })

setGeneric(name="dat.files.get.one.channel",
           def=function(df,channel.no,first.sample,last.sample)
           {standardGeneric("dat.files.get.one.channel")}
)
setMethod(f="dat.files.get.one.channel",
          signature="DatFiles",
          definition=function(df,channel.no,first.sample,last.sample)
          {
            if(length(df@file.names)==0)
              stop("df@file.names length = 0")
            if(df@num.channels<=0)
              stop("df@num.channles should be larger than 0")
            if(channel.no<0)
              stop("channel.no < 0")
            if(channel.no>=df@num.channels)
              stop("channel.no >= df@num.channels")
            if(first.sample<0)
              stop("first.sample<0")
            if(last.sample<0)
              stop("last.sample<0")
            if(first.sample>last.sample)
              stop("first.sample > last.sample")
            if(length(df@file.names)>100)
              stop("can only works with a maximum of 100 files")
            if(any(nchar(df@file.names)>255))
              stop("max file name size is 255")
            
            results<-.Call("group_data_file_si_get_one_channel_cwrap",
                df@file.names,
                df@num.channels,
                channel.no,
                first.sample,
                last.sample)
          return(results)
          })


### show ###
setMethod("show", "DatFiles",
          function(object){
            print(paste("file.names:"))
            print(paste(object@file.names))
            print(paste("num.channels:",object@num.channels))
          })



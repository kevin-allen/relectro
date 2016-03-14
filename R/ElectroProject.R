############################################
### definition of ElectroProject Class   ###
############################################
ElectroProject <- setClass(
  "ElectroProject", ## name of the class
  slots=c(directory="character",
          nSessions="numeric",
          sessionNameList="character",
          sessionList="list"),
  prototype = list(directory=""))

### setSessionList ###
setGeneric(name="setSessionList",
           def=function(ep)
           {standardGeneric("setSessionList")}
)
setMethod(f="setSessionList",
          signature="ElectroProject",
          definition=function(ep)
          {
            if(ep@directory=="")
              stop("ep@directory not set")
            
            ## list all directories in the project path
            dirs<-list.dirs(path=ep@directory)
            
            ## only keep the directory with a hyphen in the name
            dirs<-dirs[grepl(pattern="-",dirs)]
            ep@nSessions<-length(dirs)
            
            dirDepth<-length(strsplit(dirs,split="/")[[1]])
            ep@sessionNameList<-unlist(strsplit(dirs,split="/"))[seq(from=dirDepth,to=dirDepth*ep@nSessions,by=dirDepth)]
              
            for (i in 1:length(dirs))
            {
              ep@sessionList[i]<-new("RecSession",session=ep@sessionNameList[i],path=dirs[i])
            }
          ep@sessionList<-lapply(ep@sessionList,loadRecSession)  
          return(ep)
          })


### show ###
setMethod("show", "ElectroProject",
          function(object){
            print(paste("directory:",object@directory))
            print(paste("nSessions:",object@nSessions))
            print(paste("sessionNameList:"))
            print(object@sessionNameList)
            print(paste("Clustered sessions:",sum(sapply(ep@sessionList,getIsClustered))))
          })


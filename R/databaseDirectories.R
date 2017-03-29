#' Create the directories of a recording session.
#' 
#' The subject name is obtain from the session name.
#' A directory for the subject is created if it does not exist. 
#' This subject folder will be located in the directory indicated by the path argument. 
#' A directory for the session is created in the subject directory.
#' 
#' All the data (binary and text files) from a recording session are saved in the same session folder.
#' 
#' If the linkTo is set to a valid path (e.g. /data/processing), a symbolic link from the 
#' subject directory to the value of linkTo will be created.
#' 
#' By default, all permissions will be given to the directories ('777')
#' 
#' The session names should be in the format animalName-date-firstTrialLastTrial.
#' For example ka1111-31012017-0103
#' 
#' @param session Character vector giving the session name.
#' @param path Path of the root directory of the database. For example /ext_drives/d52/data/processing
#' @param linkTo A location where to put a symbolic link to the subject directory. For example /data/processing
databaseDirectories<-function(session,
                              path,
                              linkTo=NA)
  {
  if(session=="")
    stop(paste("databaseDirectories, session == \"\""))
  if(path=="")
    stop(paste("databaseDirectories, path = \"\""))
  
  splitSession<-unlist(strsplit(session,"-"))
  if(length(splitSession)!=3)
    stop(paste("databaseDirectories, the session name should contain 2 hyphens"))
  if(!dir.exists(path))
    stop(paste("databaseDirectories, path (",path,"does not exists"))  
  
  subject=splitSession[1]
  subjectDir=paste(path,subject,sep="/")
  sessionDir=paste(path,subject,session,sep="/")
  
  ## create the subject directory if not there
  if(!dir.exists(subjectDir)){
    print(paste("create",subjectDir))
    if(dir.create(subjectDir,showWarnings = TRUE)==FALSE)
      stop(paste("databaseDirectories, problem creating",subjectDir))
    Sys.chmod(subjectDir,mode='775',use_umask = FALSE) # give all permissions on the directory
  }
  ## create the symbolic link to subject directory
  if(!is.na(linkTo)){
    linkToFile<-paste(linkTo,subject,sep="/")
    if(!dir.exists(linkTo))
      stop(paste("databaseDirectories,",linkTo,"does not exist"))
    if(!file.exists(linkToFile)){
      print(paste("create symbolic link from ",subjectDir,"to",linkToFile))
      if(file.symlink(from=subjectDir,to=linkToFile)==F)
        stop(paste("databaseDirectories, problem creating",linkToFile))
      Sys.chmod(linkToFile,mode='775',use_umask = FALSE) # give all permissions on the link
    }else{
      print(paste("databaseDirectories, file",linkToFile,"already exists"))
    }
  }
  
  ## create the session directory
  if(file.exists(sessionDir)){
    stop(paste("databaseDirectories,",sessionDir,"already exists"))
  }
  print(paste("create",sessionDir))
  if(dir.create(sessionDir,showWarnings = TRUE)==FALSE)
    stop(paste("databaseDirectories, problem creating",sessionDir))
  Sys.chmod(sessionDir,mode='775',use_umask = FALSE) # give all permissions on the directory
}

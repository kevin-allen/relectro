####################################################
### functions to work with KiloSort ################
####################################################

#' Create the configurations files to run KiloSort
#' 
#' @param rs RecSession object
writeKiloSortConfigurationFiles<-function(rs,
                                          KiloSortPath="~/repo/KiloSort",
                                          npyMatlabPath="~/repo/npy-matlab")
{
  if(!dir.exists(rs@path))
     stop("rs@path does not exist")
  KiloSortPath="~/repo/KiloSort"
  npyMatlabPath="~/repo/npy-matlab"
  if(!dir.exists(KiloSortPath))
    stop(paste(KiloSortPath,"does not exist"))
  if(!dir.exists(npyMatlabPath))
    stop(paste(npyMatlabPath,"does not exist"))
  
  #########################
  ## masterFileKiloSort ###
  #########################
  print(paste("create",paste(rs@path,"masterFileKiloSort.m",sep = "/")))
  fileConn<-file(paste(rs@path,"masterFileKiloSort.m",sep = "/"))
  writeLines(c(paste("addpath(genpath('",KiloSortPath, "'))",sep=""),
             paste("addpath(genpath('",npyMatlabPath, "'))",sep=""),
             "configKiloSort;",
             "createChannelMapFile;",
             "tic;",
             "if ops.GPU",
             "    gpuDevice(1);",
             "end",
             "[rez, DATA, uproj] = preprocessData(ops);",
             "rez                = fitTemplates(rez, DATA, uproj);",
             "rez                = fullMPMU(rez, DATA);",
             "%     rez = merge_posthoc2(rez);",
             "save(fullfile(ops.root,  'rez.mat'), 'rez', '-v7.3');",
             "rezToPhy(rez, ops.root);",
             "delete(ops.fproc);"),fileConn)
  close(fileConn)

  print(paste("create",paste(rs@path,"configKiloSort.m",sep = "/")))
  fileConn<-file(paste(rs@path,"configKiloSort.m",sep = "/"))
  writeLines(c("ops.GPU = 1;",
               "ops.parfor= 0;",
               "ops.verbose= 1;",
               "ops.showfigures= 1;",
               "ops.datatype='dat';",
               paste("ops.fbinary=",paste("'",rs@fileBase,".dat';",sep="")), # one dat file per session!
               paste("ops.fproc=",paste("'",rs@fileBase,".temp_wh.dat';",sep="")),
               paste("ops.root=",paste("'",rs@path,"';",sep="")),
               paste("ops.fs=",rs@samplingRate,";"),
               paste("ops.NchanTOT=",rs@nChannels,";"),
               paste("ops.Nchan=",sum(!is.na(rs@channelsTetrode)),";"),
               "ops.Nfilt=64;",
               "ops.nNeighPC=12;",
               "ops.nNeigh=16;",
               "ops.whitening='full';",
               "ops.nSkipCov=1;",
               "ops.whiteningRange=32;",
               paste("ops.chanMap=",paste(sep="","'",paste(rs@path,"chanMap.mat",sep="/"),"';")),
               "ops.criterionNoiseChannels=0.2;",
               "ops.Nrank=3;",
               "ops.nfullpasses=6;",
               "ops.maxFR=20000;",
               "ops.fshigh=300;",
               "%ops.fslow=2000;",
               "ops.ntbuff=64;",
               "ops.scaleproc=200;",
               "ops.NT=32*1024+ops.ntbuff;", # adjust to memory 
               "ops.Th=[4 10 10];",
               "ops.lam=[5 20 20];",
               "ops.nannealpasses=4;",
               "ops.momentum=1./[20 400];",
               "ops.shuffle_clusters=1;",
               "ops.mergeT=0.1;",
               "ops.splitT=0.1;",
               "ops.initialize='no';",
               "ops.spkTh=-5;",
               "ops.loc_range=[3 1];",
               "ops.long_range=[30 6];",
               "ops.maskMaxChannels=5;",
               "ops.crit=0.65;",
               "osp.nFiltMax=10000;",
               "dd = load('PCspikes2.mat');",
               "ops.wPCA=dd.Wi(:,1:7);",
               "ops.fracse = 0.1;",
               "ops.epu=Inf;",
               "ops.ForceMaxRAMforDat=20e9;"),fileConn)
  close(fileConn)
               
  ############################
  ### createChannelMapFile ###
  ############################
  
  ## get the tetrode no of each channel
  x<-getTetrodeNumberOfChannel(rs,0:(rs@nChannels-1))
  x[which(is.na(x))]<-NaN
  nonConnected=which(!(0:(rs@nChannels-1) %in% rs@channelsTetrode))
  print(paste("create",paste(rs@path,"createChannelMapFile.m",sep = "/")))
  fileConn<-file(paste(rs@path,"createChannelMapFile.m",sep = "/"))
  writeLines(c(paste("Nchannels =",rs@nChannels,";"),
               paste("chanMap= [", paste(1:rs@nChannels,collapse=" "),"];"),
               paste("chanMap0ind = chanMap - 1;"),
               paste("connected= true(",rs@nChannels,",1);"),
               paste("connected([",paste(nonConnected,collapse=" "),"]) = 0;"),
               paste("xcoords = 20 * [", paste(x,collapse=" "),"];"),
               paste("ycoords = 20 * [", paste(rep(1,rs@nChannels),collapse = " "),"];"),
               paste("kcoords = [", paste(x,collapse=" "),"];"),
               paste("fs=",rs@samplingRate,";",sep=""),
               paste("save('", paste(rs@path,"chanMap.mat",sep="/") ,"', 'chanMap','connected','xcoords','ycoords','kcoords','chanMap0ind','fs');",sep="")),
             fileConn)
  close(fileConn)
}
  

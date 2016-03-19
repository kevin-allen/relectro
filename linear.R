## find the data provided with the package as example
clufile<-unlist(strsplit(x=system.file("extdata", "jp4298-15022016-0106.clu", package = "relectro"), split="/"))
datadir<-paste(clufile[1:length(clufile)-1],sep="/",collapse = "/")
session=strsplit(clufile[length(clufile)],split="\\.")[[1]][1]


st<-new("SpikeTrain",session=session,path=datadir)
st<-loadSpikeTrain(st)
pt<-new("Positrack",session=session,path=datadir)
pt<-loadPositrack(pt)
rs<-new("RecSession",session=session,path=datadir)
rs<-loadRecSession(rs)

# linear track data
ptlt<-setInvalidOutsideInterval(pt,s=getIntervalsEnvironment(rs,env="lt"))
ptlt<-linearzeLinearTrack(ptlt)

## make linear place maps
sp1<-new("SpatialProperties1d",session=session)

###
### create a skeleton for the package
###
###package.skeleton(name="relectro", code_files=c("~/repo/r_packages/my_src/SpikeTrain.R","~/repo/r_packages/my_src/JoinIntervals.R"),path=c("~/repo/r_packages/"))
### copy the .c and makefile files in a src directory 
###
### NAMESPACE
##
## to install from source, go in the directory where the relectro directory is located
## R CMD build relectro
## R CMD INSTALL relectro
##

library(rbenchmark) # useful to test speed, or system.time
#library("devtools")


### choose one of the two options
## option 1
## if loading the package
library(relectro)
# option 2
## if working directly from source files, load all you need, dyn.load to get the c files.
source("~/repo/r_packages/relectro/R/SpikeTrain.R")
source("~/repo/r_packages/relectro/R/Math.R")
source("~/repo/r_packages/relectro/R/JoinIntervals.R")
source("~/repo/r_packages/relectro/R/RecSession.R")
source("~/repo/r_packages/relectro/R/Positrack.R")
source("~/repo/r_packages/relectro/R/SpatialProperties2d.R")
source("~/repo/r_packages/relectro/R/DatFiles.R")
dyn.load("~/repo/r_packages/relectro/src/relectro.so")


#########################################
#### EXAMPLES RANDOM DATA SpikeTrain ####
########################################
## generate spikes for 3 neurons  
res1<-cumsum(rpois(n=100,lambda=10))
res2<-cumsum(rpois(n=200,lambda=15))
res3<-cumsum(rpois(n=300,lambda=10))
clu<-c(rep(1,100),rep(2,200),rep(3,300))
df<-data.frame(res=c(res1,res2,res3),clu=clu)
df<-df[order(df$res),] # sort according to res values
## create a SpikeTrain object from random spikes ###
st<-new("SpikeTrain")
## set the spike trains in the object
st<-setSpikeTrain(st=st,res=df$res,clu=df$clu,sampling.rate=20000)
## print object
st
## get the spike-time autocorrelation
auto<-spikeTimeAutocorrelation(st,bin.size.ms=1,window.size.ms=200)
## plot the autocorrelation
plot(auto$count,type='l',ylab="Spike count", xlab="time (ms)")
## get the mean firing rate
meanFiringRate(st)
## set time intervals to work on
st<-setIntervals(st,s=c(0,2000),e=c(1000,3000))
st
## get the mean firing rate within these intervals
meanFiringRate(st)
## set some events, spikes of clu 2
st<-setEvents(st,events=st@res[which(st@clu==2)])
cc<-spikeTimeCrosscorrelationEvents(st)
plot(cc$count[which(cc$clu==2)],type='l')
rm(clu,res1,res2,res3,df,auto,st,cc)




######################################
## test spikeTimeCrosscorrelation ####
## generate spikes for 1 neurons  ####
######################################
res1<-cumsum(rpois(n=10000,lambda=4))
clu<-c(rep(1,10000))
df<-data.frame(res=c(res1),clu=clu)
df<-df[order(df$res),] # sort according to res values
## create a SpikeTrain object from random spikes ###
st<-new("SpikeTrain")
## set the spike trains in the object
st<-setSpikeTrain(st=st,res=df$res,clu=df$clu,sampling.rate=20000)
## set some events, spikes of clu 2
st<-setEvents(st,events=st@res[which(st@clu==1)])
cc<-spikeTimeCrosscorrelationEvents(st,bin.size.ms=0.05 ,window.size.ms = 40)
plot(cc$time,cc$count,type='l',ylim=c(0,max(cc$count)))
rm(res1,clu,df,st,cc)


##############################
#### EXAMPLE JOIN INTERVAL ###
##############################
s1<-c(0,30,50,100)
e1<-c(10,46,55,150)
s2<-c(45,125)
e2<-c(55,155)
join.intervals(s1,e1,s2,e2)
rm(s1,s2,e1,e2)

#######################################
#### EXAMPLE SPIKE TRAINS FROM FILE ###
#######################################
setwd("~/repo/r_packages/data")
session="jp4298-15022016-0106"
st<-new("SpikeTrain") # create SpikeTrain object
st@session=session # set session name
st<-loadSpikeTrain(st) # load res clu and sampling rate
cross<-spikeTimeCrosscorrelation(st,bin.size.ms=1,window.size.ms = 200,probability = T) ## calculate spike-time crosscorrelation
plot(cross$prob[which(cross$clu1==2&cross$clu2==5)],ylim=c(0,max(cross$prob[which(cross$clu1==2&cross$clu2==5)])),type='l') ## plot one crosscorrelation
## set some events, in this case the spikes of clu 2
st<-setEvents(st,events=st@res[which(st@clu==2)])
cc<-spikeTimeCrosscorrelationEvents(st)
plot(cc$count[which(cc$clu==6)],ylim=c(0,max(cc$count[which(cc$clu==6)])),type='l')
rm(session,st,cross,cc)


################################
#### EXAMPLE OF MAKING PAIRS ###
################################
make.pairs(1:5)


#############################
#### EXAMPLE RecSession #####
#############################
setwd("~/repo/r_packages/data")
session="jp4298-15022016-0106"
rs<-new("RecSession")
rs@session<-session
## load configuration data from file
rs<-loadRecSession(rs)
## print session
rs
## get intervals with a given environment
getIntervalsEnvironment(rs,"lt")

#######################################
## Using SpikeTrain with RecSession ###
#######################################
setwd("~/repo/r_packages/data")
session="jp4298-15022016-0106"
st<-new("SpikeTrain",session=session)
st<-loadSpikeTrain(st)
rs<-new("RecSession",session=session)
rs<-loadRecSession(rs)
# set the intervals in spike train using output matrix from recSession
st<-setIntervals(st,s=getIntervalsEnvironment(rs,env="sqr70"))
# print the information of the new spike train
st
rm(session,st,rs)


###########################
#### Positrack example ####
###########################
setwd("~/repo/r_packages/data")
session="jp4298-15022016-0106"
pt<-new("Positrack",session=session)
pt<-loadPositrack(pt)
m<-get.intervals.at.speed(pt,0,3)
m[,1:10]
rm(session,pt,m)

######################################
### use RecSession with Positrack  ###
######################################
setwd("~/repo/r_packages/data")
session="jp4298-15022016-0106"
pt<-new("Positrack",session=session)
pt<-loadPositrack(pt)
rs<-new("RecSession",session=session)
rs<-loadRecSession(rs)
pt1<-set.invalid.outside.interval(pt,s=getIntervalsEnvironment(rs,env="sqr70"))
plot(pt1@x,pt1@y,xlab="x (cm)",ylab="y (cm)")
rm(session,pt,rs,pt1)

####################################
### get speed at some res values ###
####################################
setwd("~/repo/r_packages/data")
session="jp4298-15022016-0106"
pt<-new("Positrack",session=session)
pt<-loadPositrack(pt)
plot(get.speed.at.res.values(pt,res=seq(100000,300000,20)),type='l')
rm(pt)


################################
### SpatialProperties object ###
################################

setwd("~/repo/r_packages/data")
session="jp4298-15022016-0106"
rs<-new("RecSession",session=session) ## info about rec session
rs<-loadRecSession(rs)
pt<-new("Positrack",session=session) ## info about position
pt<-loadPositrack(pt)
st<-new("SpikeTrain",session=session) ## info about spike trains
st<-loadSpikeTrain(st)
## now do some spatial analysis with spike trains and positrack data
sp<-new("SpatialProperties2d",session=session) ## object to get spatial properties
pt<-set.invalid.outside.interval(pt,s=getIntervalsEnvironment(rs,env="sqr70")) ## select position data for one environment
sp<-firing.rate.map.2d(sp,st,pt) ## make firing rate maps
sp<-get.map.stats(sp) ## get info score, sparsity from maps
sp<-map.spatial.autocorrelation(sp) ## spatial autocorrelation from maps
sp<-grid.score(sp)
sp<-grid.orientation(sp)
sp<-grid.spacing(sp)
sp<-border.score(sp)
## plot one map
jet.colors = colorRampPalette(c("#00007F", "blue","#007FFF",  "cyan", "#7FFF7F", "yellow", "#FF7F00","red"))
map<-sp@maps[,,7]
image(t(map),zlim=c(0,max(map,na.rm=T)), col=jet.colors(200),xlab='',ylab='',axes=FALSE)
rm(sp,pt,jet.colors,map,st,session,rs)



########################
### DatFiles object ####
########################
## get session information
rs<-new("RecSession")
setwd("/data/processing/jp4298/jp4298-12022016-0104/")
rs@session="jp4298-12022016-0104"
rs<-loadRecSession(rs)
## get data from files, use rs to get info we need
df<-new("DatFiles")
df<-dat.files.set(df,file.names=paste(rs@trial.names,"dat",sep="."),num.channels=rs@n.channels)
data<-dat.files.get.one.channel(df,channel.no=df@num.channels-1,first.sample=0,last.sample=rs@trial.end.res[length(rs@trial.end.res)]) # get data one channel
## detect ttl ups or down in signal
ups<-detect.ups(x=data)
downs<-detect.downs(x=data)
plot(ups[1:50],downs[1:50])
lines(c(0,120000),c(0,120000))
rm(rs,df,downs,ups,session,data)

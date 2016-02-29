###
### create a skeleton for the package
###
###package.skeleton(name="relectro", code_files=c("~/repo/r_packages/my_src/SpikeTrain.R","~/repo/r_packages/my_src/JoinIntervals.R"),path=c("~/repo/r_packages/"))
## copy the .c and makefile files in a src directory 
##
## NAMESPACE

##library(rbenchmark)
#library("devtools")

## load the package
library(relectro)


## to reload a modified package, restart R session and call library again
#load_all("relectro",recompile=T)

## when working directly from source files, load all you need, dyn.load to get the c files.
source("~/repo/r_packages/relectro/R/SpikeTrain.R")
source("~/repo/r_packages/relectro/R/Math.R")
source("~/repo/r_packages/relectro/R/JoinIntervals.R")
source("~/repo/r_packages/relectro/R/RecSession.R")
dyn.load("~/repo/r_packages/relectro/src/relectro.so")


#########################################
#### EXAMPLES RANDOM DATA SpikeTrain ####
########################################
## generate spikes for 2 neurons  
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
plot(auto$count,type='l')
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

######################################################################
## set intervals of SpikeTrain for a a given recording environment ###
######################################################################
setwd("~/repo/r_packages/data")
session="jp4298-15022016-0106"
# get the SpikeTime
st<-new("SpikeTrain")
st@session=session
st<-loadSpikeTrain(st)
# get the RecSession
rs<-new("RecSession")
rs@session<-session
rs<-loadRecSession(rs)
# set the intervals in spike train using output matrix from recSession
st<-setIntervals(st,s=getIntervalsEnvironment(rs,env="lt"))
# print the information of the new spike train
rs
rm(session,st,rs)


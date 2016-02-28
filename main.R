###
### create a skeleton for the package
###
###package.skeleton(name="relectro", code_files=c("~/repo/r_packages/my_src/SpikeTrain.R","~/repo/r_packages/my_src/JoinIntervals.R"),path=c("~/repo/r_packages/"))
## copy the .c and makefile files in a src directory 
##
## NAMESPACE

##library(rbenchmark)
library("devtools")

## load the package
library(relectro)

## to reload a modified package, restart R session and call library again
load_all("relectro",recompile=T)

## when working directly from source files, use dyn.load to get the c files.
##dyn.load("~/repo/r_packages/relectro/src/relectro.so")
##dyn.load("~/repo/r_packages/relectro/src/firing_rate.so")


##############################
#### EXAMPLES RANDOM DATA ####
##############################
## generate spikes for 2 neurons  
res1<-cumsum(rpois(n=100,lambda=10))
res2<-cumsum(rpois(n=200,lambda=15))
clu<-c(rep(1,100),rep(2,200))
df<-data.frame(res=c(res1,res2),clu=clu)
df<-df[order(df$res),] # sort according to res values
## create a SpikeTrain object from random spikes ###
st<-new("SpikeTrain")
## set the spike trains in the object
st<-setSpikeTrain(st=st,res=df$res,clu=df$clu,sampling.rate=20000)
## get the spike-time autocorrelation
auto<-spikeTimeAutocorrelation(st,bin.size.ms=1,window.size.ms=200)
## get the mean firing rate
meanFiringRate(st)
## plot the autocorrelation
plot(auto$count,type='l')
rm(clu,res1,res2,df,auto,st)

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
rm(session,st)


################################
#### EXAMPLE OF MAKING PAIRS ###
################################
make.pairs(1:2)




##################################################################################
##### WKPETSAMP2 extension of WGBYC 2022 simulation code
##### 7 Mar 2023
##### David Lusseau
##################################################################################
###WGBYC code


library(extraDistr)
### a fishing year

fishing.day<-1:365

#fleet
nboat<-100
nmetier<-1
fleet<-1:nboat
metier<-1:nmetier; nyears <- 20 #used to be 100

#n.fishing.event<-100000


large.events<-c(0,0.0001,0.001,.005,0.01,0.05,.1,1)
p_monitor<-c(0.005,.01,0.05,.1)


p.bycatch<-.1

mean.bycatch.event<-1
mean.bycatch.large.event<-20
p.large.event<-0.01

#simple case the number 
event.type<-rbinom(1,1,p.large.event)

max((1-event.type)*rtpois(10,mean.bycatch.event,a=0),event.type*rtpois(10,mean.bycatch.large.event,a=0))

mean.fishing.event.boat.day<-2

#grow the fishing year - will optimise this process 

#initialise
#fishing.day 1


make_fishing_year<-function(mean.bycatch.event=1,mean.bycatch.large.event=20,p.large.event=0.01,
							nboat=100,mean.fishing.event.boat.day=2,p.bycatch=0.1,stochastic=TRUE) {
#we first only deal with one metier at a time
#bycatch is not affected by vessel characteeristics
#later we can for example introduce vessel size for each boat 
# and influence probabilities by vessel size
#and subsequently influence monitoring by vessel size

fishing.day<-1:365
fleet<-1:nboat							
if (stochastic==TRUE) {

mean.fishing.event.boat.day<-rtpois(nboat,mean.fishing.event.boat.day,a=0)  #introduce stochasticity so that the mean number of events per boats vary
fishing.event.per.boat<-rpois(nboat,mean.fishing.event.boat.day)

} else {
fishing.event.per.boat<-rpois(nboat,mean.fishing.event.boat.day) #uniform fishing behaviour

}

i=1
fishing<-data.frame(fishing.day=fishing.day[i],boat=rep(fleet,fishing.event.per.boat),bycatch=rbinom(sum(fishing.event.per.boat),1,p.bycatch),nbycatch=0)
event.type<-rbinom(sum(fishing$bycatch),1,p.large.event)
fishing$nbycatch[fishing$bycatch==1]<-apply(cbind((1-event.type)*rtpois(sum(fishing$bycatch),mean.bycatch.event,a=0),event.type*rtpois(sum(fishing$bycatch),mean.bycatch.large.event,a=0)),1,max)


for (i in 2:365) {

if (stochastic==TRUE) {
mean.fishing.event.boat.day<-rtpois(nboat,mean.fishing.event.boat.day,a=0)  #introduce stochasticity so that the mean number of events per boats vary
fishing.event.per.boat<-rpois(nboat,mean.fishing.event.boat.day)
} else {
fishing.event.per.boat<-rpois(nboat,mean.fishing.event.boat.day) #uniform fishing behaviour
}

temp<-data.frame(fishing.day=fishing.day[i],boat=rep(fleet,fishing.event.per.boat),bycatch=rbinom(sum(fishing.event.per.boat),1,p.bycatch),nbycatch=0)
event.type<-rbinom(sum(temp$bycatch),1,p.large.event)
temp$nbycatch[temp$bycatch==1]<-apply(cbind((1-event.type)*rtpois(sum(temp$bycatch),mean.bycatch.event,a=0),event.type*rtpois(sum(temp$bycatch),mean.bycatch.large.event,a=0)),1,max)

fishing<-rbind(fishing,temp)

}

return(fishing)
}


#first we simply look at the effect of monitoring rate
##################################################################################
###################################################################################
#### functions

monitor_BPUE<-function(pmonitor=0.5,nsample=1000,BPUE_real=0,fishing=NA, p_monitor_boat=.1,boat_samp=TRUE) {

#integrate observer procedure effect
#this is protocol influence process
# two potential effects : haul level effect: unobserved bycatch when monitored (lack of observation)
#as well as decrease in the number of animal reported on a bycatch event (incomplete observation, eg drop out, multiple location needed)
#observation probability function post sampling

#second step: variance associated with observed ID


BPUE_est<-array(nsample)

for (i in 1:nsample) {

if (boat_samp==TRUE) {
monitored<-sample(c(1:dim(fishing)[1]),floor(pmonitor*dim(fishing)[1]),replace=FALSE) # sample without replacement
BPUE_est[i]<-(sum(fishing$nbycatch[monitored])/length(monitored))
} else {
boat_sampled<-sample(1:max(fishing$boat),n=floor(max(fishing$boat)*p_monitor_boat),replace=FALSE)
fleet_sampled<-fishing[fishing$boat%in%boat_sampled,]
monitored<-sample(c(1:dim(fleet_sampled)[1]),floor(pmonitor*dim(fleet_sampled)[1]),replace=FALSE) # sample without replacement
BPUE_est[i]<-(sum(fleet_sampled$nbycatch[monitored])/length(monitored))
}

}
BPUE_est_mean<-mean(BPUE_est)
BPUE_est_CV<-sd(BPUE_est)/mean(BPUE_est)

return(list(BPUE_est=BPUE_est_mean,CV=BPUE_est_CV))
}


#extension:
#only monitor a small proportion of the fleet




fishing1<-make_fishing_year()

p_monitor<-c(seq(.01,.2,.01),seq(.25,.5,.05))

BPUE_real<-sum(fishing1$nbycatch)/dim(fishing1)[1]

#try 1 monitoring programme
BPUE_estimate<-monitor_BPUE(pmonitor=p_monitor[1],BPUE_real=BPUE_real,fishing=fishing1)

## first assessment what happens when I increase monitoring rate

monitor_estimate<-data.frame(p_monitor=p_monitor,BPUE_real=BPUE_real,BPUE_est=NA,BPUE_est_CV=NA)

for (i in 1:length(p_monitor)) {

monitor<-monitor_BPUE(pmonitor=p_monitor[i],BPUE_real=BPUE_real,fishing=fishing1)
monitor_estimate$BPUE_est[i]<-monitor$BPUE_est
monitor_estimate$BPUE_est_CV[i]<-monitor$CV
}

## now we are going to replicate the real fishing year replication 100 times
monitor_estimate$rel.bias<-(monitor_estimate$BPUE_real-monitor_estimate$BPUE_est)/monitor_estimate$BPUE_real
monitor_estimate$year<-1

for (y in 1:nyears) {
fishing1<-make_fishing_year()
BPUE_real<-sum(fishing1$nbycatch)/dim(fishing1)[1]


temp<-data.frame(p_monitor=p_monitor,BPUE_real=BPUE_real,BPUE_est=NA,BPUE_est_CV=NA)

for (i in 1:length(p_monitor)) {

monitor<-monitor_BPUE(pmonitor=p_monitor[i],BPUE_real=BPUE_real,fishing=fishing1)
temp$BPUE_est[i]<-monitor$BPUE_est
temp$BPUE_est_CV[i]<-monitor$CV
}

temp$rel.bias<-(temp$BPUE_real-temp$BPUE_est)/temp$BPUE_real
temp$year<-y


monitor_estimate<-rbind(monitor_estimate,temp)

}


##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################

##############################################################################################################################
#### 7 March 2023
#### first premise, effort reporting is not complete
#### what happens first if we are randomly not reporting fishing events from some vessels and/or from some hauls
#### the next extension on this will make it inaccurate



library(extraDistr)
### a fishing year

fishing.day<-1:365

#fleet
nboat<-100
nmetier<-1
fleet<-1:nboat
metier<-1:nmetier

#n.fishing.event<-100000


large.events<-c(0,0.0001,0.001,.005,0.01,0.05,.1,1)
p_monitor<-c(0.005,.01,0.05,.1)


p.bycatch<-.1

mean.bycatch.event<-1
mean.bycatch.large.event<-20
p.large.event<-0.01

#simple case the number 
event.type<-rbinom(1,1,p.large.event)

max((1-event.type)*rtpois(10,mean.bycatch.event,a=0),event.type*rtpois(10,mean.bycatch.large.event,a=0))

mean.fishing.event.boat.day<-2

#grow the fishing year - will optimise this process 

#initialise
#fishing.day 1


make_fishing_year<-function(mean.bycatch.event=1,mean.bycatch.large.event=20,p.large.event=0.01,
							nboat=100,mean.fishing.event.boat.day=2,p.bycatch=0.1,stochastic=TRUE) {
#we first only deal with one metier at a time
#bycatch is not affected by vessel characteeristics
#later we can for example introduce vessel size for each boat 
# and influence probabilities by vessel size
#and subsequently influence monitoring by vessel size

fishing.day<-1:365
fleet<-1:nboat							
if (stochastic==TRUE) {

mean.fishing.event.boat.day<-rtpois(nboat,mean.fishing.event.boat.day,a=0)  #introduce stochasticity so that the mean number of events per boats vary
fishing.event.per.boat<-rpois(nboat,mean.fishing.event.boat.day)

} else {
fishing.event.per.boat<-rpois(nboat,mean.fishing.event.boat.day) #uniform fishing behaviour

}

i=1
fishing<-data.frame(fishing.day=fishing.day[i],boat=rep(fleet,fishing.event.per.boat),bycatch=rbinom(sum(fishing.event.per.boat),1,p.bycatch),nbycatch=0)
event.type<-rbinom(sum(fishing$bycatch),1,p.large.event)
fishing$nbycatch[fishing$bycatch==1]<-apply(cbind((1-event.type)*rtpois(sum(fishing$bycatch),mean.bycatch.event,a=0),event.type*rtpois(sum(fishing$bycatch),mean.bycatch.large.event,a=0)),1,max)



#########
## so for this challenge we need to change the computation of the estimated total bycatch it becomes the estimated BPUE x estimated effort
return(fishing)
}

estimate_fishing_effort<-function(fishing=NA,p_report=.9) {
effort<-data.frame(DaS_report=NA,hauls_report=NA,DaS_estimated=NA,hauls_estimated=NA,DaS_real=NA,hauls_real=NA)
#we calculate real effort, reported effort, and estimated effort where we know the fleet size and assume average effort for those unreported
# here reporting bias is at the vessel/boat level at first random then biased by high res metier (eg small vessels not reporting)
reporting<-sample(1:nboat,floor(p_report*nboat),replace=FALSE)
fishing.report<-fishing[fishing$boat%in%reporting,]
effort$hauls_real<-dim(fishing)[1]
effort$hauls_report<-dim(fishing.report)[1]
effort$hauls_estimated<-dim(fishing.report)[1]+(nboat-floor(p_report*nboat))*mean(table(fishing.report$boat)) #add 1-p_report nboat fishing mean number of hauls of those reporting
Das<-rowSums(table(fishing$boat,fishing$fishing.day)!=0)
Das_report<-rowSums(table(fishing.report$boat,fishing.report$fishing.day)!=0)

effort$DaS_real<-sum(Das)
effort$DaS_report<-sum(Das_report)
effort$DaS_estimated=sum(Das_report)+(nboat-floor(p_report*nboat))*mean(Das)

return(effort)
}

#################################################################################################################
#################################################################################################################
# first look at effort
p_reported<-seq(0.1,1,.1)

effort_bias<-data.frame(year=rep(1:100,each=length(p_reported)),reporting=rep(p_reported,100),DaS_report=NA,hauls_report=NA,DaS_estimated=NA,hauls_estimated=NA,DaS_real=NA,hauls_real=NA)

m<-1
for (j in 1:100) { #replicate 100 years
fishing<-make_fishing_year()

for (i in 1:length(p_reported)) {

effort_bias[m,3:8]<-estimate_fishing_effort(fishing=fishing,p_report=p_reported[i])

m<-m+1
print(m)
flush.console()
}
}

library(ggplot2)
effort_bias$DaS_report_bias<-(effort_bias$DaS_real-effort_bias$DaS_report)/effort_bias$DaS_real
effort_bias$DaS_estimated_bias<-(effort_bias$DaS_real-effort_bias$DaS_estimated)/effort_bias$DaS_real
effort_bias$hauls_report_bias<-(effort_bias$hauls_real-effort_bias$hauls_report)/effort_bias$hauls_real
effort_bias$hauls_estimated_bias<-(effort_bias$hauls_real-effort_bias$hauls_estimated)/effort_bias$hauls_real

ggplot(effort_bias)+
geom_point(aes(x=reporting,y=DaS_report_bias),colour="blue")+
geom_point(aes(x=reporting,y=DaS_estimated_bias),colour="red")+
labs(x="proportion of the fleet reporting effort",y="Days at Sea")


ggplot(effort_bias,aes(x=reporting,y=DaS_estimated_bias))+
geom_point()+
labs(x="proportion of the fleet reporting effort",y="Days at Sea")

ggplot(effort_bias,aes(x=reporting,y=DaS_estimated))+
geom_point()+
labs(x="proportion of the fleet reporting effort",y="Days at Sea")




#################################################################################################################
#################################################################################################################
###total mortality
#BPUE(all forms) x effort (all forms)

#################################################################################################################
#################################################################################################################

#### introducing varying boat characteristics which influence bycatch probability
## treating at first as two 'high res' metier categories, so say metier level 6 for a fishing year inspected at level 4

make_fishing_year_metier<-function(mean.bycatch.event=1,mean.bycatch.large.event=20,p.large.event=0.01,
							nboat=100,mean.fishing.event.boat.day=2,p.bycatch=c(0.1,.01),p.metier=c(.2,.8),stochastic=TRUE) {

#p.metier is the proportion of vessel in the, here, length of p.bycatch metiers 
# p bycatch event alternative distribution particularly for low density species

nmetier<-length(p.bycatch)

fishing.day<-1:365
fleet<-1:nboat							
metier<-sample(1:nmetier,nboat,replace=TRUE,prob=p.metier) #here we deal with probability of metier, not proportion of metier

if (stochastic==TRUE) {
# here number of hauls is still not associated to "high res" metier
mean.fishing.event.boat.day<-rtpois(nboat,mean.fishing.event.boat.day,a=0)  #introduce stochasticity so that the mean number of events per boats vary
fishing.event.per.boat<-rpois(nboat,mean.fishing.event.boat.day)

} else {
fishing.event.per.boat<-rpois(nboat,mean.fishing.event.boat.day) #uniform fishing behaviour

}

i=1
fishing<-data.frame(fishing.day=fishing.day[i],boat=rep(fleet,fishing.event.per.boat),metiers=rep(metier,fishing.event.per.boat),bycatch=rbinom(sum(fishing.event.per.boat),1,p.bycatch[rep(metier,fishing.event.per.boat)]),nbycatch=0)

event.type<-rbinom(sum(fishing$bycatch),1,p.large.event)
fishing$nbycatch[fishing$bycatch==1]<-apply(cbind((1-event.type)*rtpois(sum(fishing$bycatch),mean.bycatch.event,a=0),event.type*rtpois(sum(fishing$bycatch),mean.bycatch.large.event,a=0)),1,max)


for (i in 2:365) {

if (stochastic==TRUE) {
#mean.fishing.event.boat.day<-rtpois(nboat,mean.fishing.event.boat.day,a=0) #that's stays the same for the whole year #introduce stochasticity so that the mean number of events per boats vary
fishing.event.per.boat<-rpois(nboat,mean.fishing.event.boat.day)
} else {
fishing.event.per.boat<-rpois(nboat,mean.fishing.event.boat.day) #uniform fishing behaviour
}

temp<-data.frame(fishing.day=fishing.day[i],boat=rep(fleet,fishing.event.per.boat),metiers=rep(metier,fishing.event.per.boat),bycatch=rbinom(sum(fishing.event.per.boat),1,p.bycatch[rep(metier,fishing.event.per.boat)]),nbycatch=0)
event.type<-rbinom(sum(temp$bycatch),1,p.large.event)
temp$nbycatch[temp$bycatch==1]<-apply(cbind((1-event.type)*rtpois(sum(temp$bycatch),mean.bycatch.event,a=0),event.type*rtpois(sum(temp$bycatch),mean.bycatch.large.event,a=0)),1,max)

fishing<-rbind(fishing,temp)

}
#########
## so for this challenge we need to change the computation of the estimated total bycatch it becomes the estimated BPUE x estimated effort
return(fishing)
}

#################################################################################################################
#################################################################################################################

#################################################################################################################

#Introducing area and time variability

## treating at first as two 'high res' metier categories, so say metier level 6 for a fishing year inspected at level 4

#parameters controlling fishery effort allocation and bycatch hot spot areas, relative to time
narea<-10 #defines how many areas you want
spatial.effort.skewness.general<-c(1,1) #alpha and beta parameters from the in beta-binomial distribution, X∼BB(n,α,β). If α=β=1, then it is a discrete uniform distribution, if α≥1 and β<1 then it is a discrete left-skewed distribution.
spatial.effort.skewness.special<-c(1.7,0.3)# for time periods with different distribution
time.periods.fishery<-32:60 #time of year you want the within area effort to shift, could be several periods (uses spatial.effort.skewness.special )
time.periods.bycatch<-32:60#if similar as time.periods.fishery
hostspot.area<-10 #define area(s) that you would like to be hotspots for bycatch
##
make_fishing_year_metier<-function(mean.bycatch.event=1,mean.bycatch.large.event=20,p.large.event=0.01,
                                   nboat=100,mean.fishing.event.boat.day=2,p.bycatch=c(0.1,.01),p.metier=c(.2,.8),
                                   narea=10,stochastic=TRUE, spatio.temporal.fishery.trend=TRUE,spatio.temporal.bycatch.trend=TRUE,
                                   spatial.effort.skewness.general=c(1,1), spatial.effort.skewness.special=c(1.7,0.3),
                                   time.periods.fishery=32:60,time.periods.bycatch=32:60,hostspot.area=10
                                   ) {
  
  #p.metier is the proportion of vessel in the, here, length of p.bycatch metiers 
  # p bycatch event alternative distribution particularly for low density species
  
  nmetier<-length(p.bycatch)
  
  fishing.day<-1:365
  fleet<-1:nboat
  metier<-sample(1:nmetier,nboat,replace=TRUE,prob=p.metier) #here we deal with probability of metier, not proportion of metier
  
  if (stochastic==TRUE) {
    # here number of hauls is still not associated to "high res" metier
    mean.fishing.event.boat.day<-rtpois(nboat,mean.fishing.event.boat.day,a=0)  #introduce stochasticity so that the mean number of events per boats vary
    fishing.event.per.boat<-rpois(nboat,mean.fishing.event.boat.day)
    
  } else {
    fishing.event.per.boat<-rpois(nboat,mean.fishing.event.boat.day) #uniform fishing behaviour
    
  }
  

fishing.area.dist<-sample(1:narea,nboat,replace=TRUE,prob=dbbinom(1:narea-1, narea-1, spatial.effort.skewness.general[1], spatial.effort.skewness.general[2]))#(runif(nboat,1,narea)) #
  
    
 
  
  i=1
  fishing<-data.frame(fishing.day=fishing.day[i],boat=rep(fleet,fishing.event.per.boat),area=rep(fishing.area.dist,fishing.event.per.boat),metiers=rep(metier,fishing.event.per.boat),bycatch=rbinom(sum(fishing.event.per.boat),1,p.bycatch[rep(metier,fishing.event.per.boat)]),nbycatch=0)
  
  event.type<-rbinom(sum(fishing$bycatch),1,p.large.event)
  fishing$nbycatch[fishing$bycatch==1]<-apply(cbind((1-event.type)*rtpois(sum(fishing$bycatch),mean.bycatch.event,a=0),event.type*rtpois(sum(fishing$bycatch),mean.bycatch.large.event,a=0)),1,max)
  
    
  
  for (i in 2:365) {
    
    if (stochastic==TRUE) {
      #mean.fishing.event.boat.day<-rtpois(nboat,mean.fishing.event.boat.day,a=0) #that's stays the same for the whole year #introduce stochasticity so that the mean number of events per boats vary
      fishing.event.per.boat<-rpois(nboat,mean.fishing.event.boat.day)
    } else {
      fishing.event.per.boat<-rpois(nboat,mean.fishing.event.boat.day) #uniform fishing behaviour
    }
  
    if (spatio.temporal.fishery.trend==TRUE){
      #enable change in fishery density in specific areas and time periods
      fishing.area.dist = if(i %in% time.periods.fishery){sample(1:narea,nboat,replace=TRUE,prob=dbbinom(1:narea-1, narea-1, spatial.effort.skewness.special[1], spatial.effort.skewness.special[2]))
      } else{
      sample(1:narea,nboat,replace=TRUE,prob=dbbinom(1:narea-1, narea-1, spatial.effort.skewness.general[1], spatial.effort.skewness.general[2]))
      }
      area=rep(fishing.area.dist,fishing.event.per.boat)
      }else{
        fishing.area.dist<-sample(1:narea,nboat,replace=TRUE,prob=dbbinom(1:narea-1, narea-1, spatial.effort.skewness.general[1], spatial.effort.skewness.general[2])) 
      area=rep(fishing.area.dist,fishing.event.per.boat)
    }     
     
    
    temp<-data.frame(fishing.day=fishing.day[i],boat=rep(fleet,fishing.event.per.boat),area=area,metiers=rep(metier,fishing.event.per.boat),bycatch=rbinom(sum(fishing.event.per.boat),1,p.bycatch[rep(metier,fishing.event.per.boat)]),nbycatch=0)
      #enable change in fishery density in specific areas and time periods
    tryCatch({  
    bycatch_high = temp[temp$fishing.day %in% time.periods.bycatch & temp$area %in%  hostspot.area,]
      bycatch_high$bycatch<-bycatch=rbinom(sum(fishing.event.per.boat),1,p.bycatch[rep(metier,fishing.event.per.boat)]*2)
      bycatch_low = temp[!temp$fishing.day %in% time.periods.bycatch & !temp$area %in%  hostspot.area,]
      bycatch_low$bycatch<-bycatch=rbinom(sum(fishing.event.per.boat),1,p.bycatch[rep(metier,fishing.event.per.boat)])
      fisherydata<-rbind(bycatch_high,bycatch_low)},
  error=function(e){})
    
    event.type<-rbinom(sum(temp$bycatch),1,p.large.event)
    temp$nbycatch[temp$bycatch==1]<-apply(cbind((1-event.type)*rtpois(sum(temp$bycatch),mean.bycatch.event,a=0),event.type*rtpois(sum(temp$bycatch),mean.bycatch.large.event,a=0)),1,max)
    
    fishing<-rbind(fishing,temp)
    
  }
  #########
  ## so for this challenge we need to change the computation of the estimated total bycatch it becomes the estimated BPUE x estimated effort
  return(fishing)
}



#################################################################################################################
#################################################################################################################
# look at effort and then introducing high res metier effort reporting bias
p_reported<-seq(0.1,1,.1)

effort_bias<-data.frame(year=rep(1:100,each=length(p_reported)),reporting=rep(p_reported,100),DaS_report=NA,hauls_report=NA,DaS_estimated=NA,hauls_estimated=NA,DaS_real=NA,hauls_real=NA)

m<-1
for (j in 1:100) { #replicate 100 years
fishing<-make_fishing_year_metier()

for (i in 1:length(p_reported)) {

effort_bias[m,3:8]<-estimate_fishing_effort(fishing=fishing,p_report=p_reported[i])

m<-m+1
print(m)
flush.console()
}
}

library(ggplot2)
effort_bias$DaS_report_bias<-(effort_bias$DaS_real-effort_bias$DaS_report)/effort_bias$DaS_real
effort_bias$DaS_estimated_bias<-(effort_bias$DaS_real-effort_bias$DaS_estimated)/effort_bias$DaS_real
effort_bias$hauls_report_bias<-(effort_bias$hauls_real-effort_bias$hauls_report)/effort_bias$hauls_real
effort_bias$hauls_estimated_bias<-(effort_bias$hauls_real-effort_bias$hauls_estimated)/effort_bias$hauls_real

ggplot(effort_bias)+
geom_point(aes(x=reporting,y=DaS_report_bias),colour="blue")+
geom_point(aes(x=reporting,y=DaS_estimated_bias),colour="red")+
labs(x="proportion of the fleet reporting effort",y="Days at Sea")


ggplot(effort_bias,aes(x=reporting,y=DaS_estimated_bias))+
geom_point()+
labs(x="proportion of the fleet reporting effort",y="Days at Sea")

ggplot(effort_bias,aes(x=reporting,y=DaS_estimated))+
geom_point()+
labs(x="proportion of the fleet reporting effort",y="Days at Sea")

ggplot(effort_bias,aes(x=reporting,y=hauls_estimated))+
geom_point()+
labs(x="proportion of the fleet reporting effort",y="number of hauls")

ggplot(effort_bias,aes(x=reporting,y=hauls_estimated_bias))+
geom_point()+
labs(x="proportion of the fleet reporting effort",y="number of hauls bias")

#note the difference between the variance of DaS and hauls. here we are in a situation where given the haul/day mean, most vessels
# will go out almost everyday, so the DaS 'extrapolation' is fairly robust, but then the haul number is not because this varies
# a little between boats





#the reporting rate is associated with vessel length
#what happens if small vessels are not reported, they represent a large proportion of the fleet (say 80%) but also
# their fishing days differ: 
# option 1: fishing days contain less hauls
# option 2: they fish in high density areas (closer to shore) so their bycatch rate is higher

#### second premise we deal varying bycatch probability with metier level 6

#### third premise we deteriorate the monitoring programme strata to metier 3

#### fourth premise we deal with bycatch complexity - multiple fisheries and multiple species 
# first PTB and GNS






#######################3
#### 8 Mar 2023
#### task 1: observer process
#### task 2: within year trends in variability in fishing effort (mean event per day) and bycatch probability
#### task 3: introduce "metier" in monitoring and effort function

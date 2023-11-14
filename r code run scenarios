##################################################################################
##### WKPETSAMP2 extension of WGBYC 2022 simulation code
##### 7 Mar 2023
##### David Lusseau, Paula Gutierrez, Henrik Pärn, Kim, Torbjörn, Margaret Siple
##################################################################################
###WGBYC code


library(extraDistr)
### a fishing year


##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################

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

monitor_BPUE_metier<-function(pmonitor=0.5,nsample=1000,BPUE_real=0,fishing=NA, p_monitor_boat=.1,boat_samp=TRUE,
p_haul_obs=1,detect_prob=1,refusal_rate=0, misclassification=0, bymetier=FALSE, p_monitor_metier=.5) {

#p_haul_obs the probability that an entire fishing events (hauls) was observed by the observer
#detect_prob is the detection probability of each individual in the bycatch event
#while at first, detect_prob is an independent draw on each individual in the haul, we can implement a decrease in detect_prob
#with nbycatch increasing (on a log scale or similar)

#misclassification is the mis-identification probability for the bycaught species
#misclassification might come in later when we have multiple species which can be confused, and therefore individuals from nbycatch_sp_i can
#be taken to add to nbycatch_sp_j

#refusal_rate: the probability that a vessel selected for monitoring refuses to engage
#refusal is replaced to introduce selective pressures
#option 1 random option 2 associated with nbycatch

#integrate observer procedure effect
#this is protocol influence process
# two potential effects : haul level effect: unobserved bycatch when monitored (lack of observation)
#as well as decrease in the number of animal reported on a bycatch event (incomplete observation, eg drop out, multiple location needed)
#observation probability function post sampling

#second step: variance associated with observed ID

#p_monitor_metier: what is the proportion of the monitoring effort which is attributed to each "metier" simulated in fishing
nmetier=length(unique(fishing$metiers))

if (bymetier==FALSE) {
BPUE_est<-array(NA,dim=nsample)
} else {
BPUE_est<-array(NA,dim=c(nsample,nmetier))
}

for (j in 1:nsample) { 

if (bymetier==FALSE) {

if (boat_samp==FALSE) {
monitored<-sample(c(1:dim(fishing)[1]),floor(pmonitor*dim(fishing)[1]),replace=FALSE) # sample without replacement
fishing_monitored<-fishing[monitored,]
not_observed<-sample(c(1:dim(fishing_monitored)[1]),floor((1-p_haul_obs)*dim(fishing_monitored)[1]),replace=FALSE)
fishing_monitored$bycatch[not_observed]<-0
fishing_monitored$nbycatch[not_observed]<-0
fishing_monitored$nbycatch<-sapply(fishing_monitored$nbycatch,function(x) rbinom(1,x,detect_prob))

BPUE_est[j]<-(sum(fishing_monitored$nbycatch)/length(monitored))

} else {
boat_sampled<-sample(unique(fishing$boat),size=floor(length(unique(fishing$boat))*p_monitor_boat),replace=FALSE)
## think about situations when boats are never going out in a year. at the moment, the observer programme allows to react and sample
## only those vessels that have been fishing at least once per year.

#refusal
boat_sampled<-sample(boat_sampled,size=floor(length(boat_sampled)*(1-refusal_rate)),replace=FALSE)

fleet_sampled<-fishing[fishing$boat%in%boat_sampled,]
monitored<-sample(c(1:dim(fleet_sampled)[1]),floor(pmonitor*dim(fleet_sampled)[1]),replace=FALSE) # sample without replacement
fishing_monitored<-fleet_sampled[monitored,]
not_observed<-sample(c(1:dim(fishing_monitored)[1]),floor((1-p_haul_obs)*dim(fishing_monitored)[1]),replace=FALSE)
fishing_monitored$bycatch[not_observed]<-0
fishing_monitored$nbycatch[not_observed]<-0
fishing_monitored$nbycatch<-sapply(fishing_monitored$nbycatch,function(x) rbinom(1,x,detect_prob))

BPUE_est[j]<-(sum(fishing_monitored$nbycatch)/length(monitored))
}

} else {
nmetier<-length(unique(fishing$metiers)) # for error catching later on

if (boat_samp==FALSE) {
nmetier<-length(unique(fishing$metiers)) # for error catching later on
monitored_by_metier<-p_monitor_metier*pmonitor

monitored<-sample(as.numeric(row.names(fishing[fishing$metiers==1,])),floor(monitored_by_metier[1]*length(row.names(fishing[fishing$metiers==1,]))),replace=FALSE) # sample without replacement
for (i in 2:nmetier) {
monitored<-c(monitored,sample(as.numeric(row.names(fishing[fishing$metiers==i,])),floor(monitored_by_metier[i]*length(row.names(fishing[fishing$metiers==i,]))),replace=FALSE)) # sample without replacement

}

fishing_monitored<-fishing[monitored,]
not_observed<-sample(c(1:dim(fishing_monitored)[1]),floor((1-p_haul_obs)*dim(fishing_monitored)[1]),replace=FALSE)
fishing_monitored$bycatch[not_observed]<-0
fishing_monitored$nbycatch[not_observed]<-0
fishing_monitored$nbycatch<-sapply(fishing_monitored$nbycatch,function(x) rbinom(1,x,detect_prob))

BPUE_est[j,]<-(tapply(fishing_monitored$nbycatch,fishing_monitored$metiers,sum)/tapply(fishing_monitored$nbycatch,fishing_monitored$metiers,length))



} else {
boat_monitored_by_metier<-p_monitor_metier*p_monitor_boat


boat_sampled<-sample(unique(fishing$boat[fishing$metiers==1]),floor(boat_monitored_by_metier[1]*length(unique(fishing$boat[fishing$metiers==1]))),replace=FALSE) # sample without replacement
for (i in 2:nmetier) {
boat_sampled<-c(boat_sampled,sample(unique(fishing$boat[fishing$metiers==i]),floor(boat_monitored_by_metier[i]*length(unique(fishing$boat[fishing$metiers==i]))),replace=FALSE)) # sample without replacement

}

#refusal
boat_sampled<-sample(boat_sampled,size=floor(length(boat_sampled)*(1-refusal_rate)),replace=FALSE)

fleet_sampled<-fishing[fishing$boat%in%boat_sampled,]
monitored<-sample(c(1:dim(fleet_sampled)[1]),floor(pmonitor*dim(fleet_sampled)[1]),replace=FALSE) # sample without replacement
fishing_monitored<-fleet_sampled[monitored,]
not_observed<-sample(c(1:dim(fishing_monitored)[1]),floor((1-p_haul_obs)*dim(fishing_monitored)[1]),replace=FALSE)
fishing_monitored$bycatch[not_observed]<-0
fishing_monitored$nbycatch[not_observed]<-0
fishing_monitored$nbycatch<-sapply(fishing_monitored$nbycatch,function(x) rbinom(1,x,detect_prob))


BPUE_est[j,]<-(tapply(fishing_monitored$nbycatch,fishing_monitored$metiers,sum)/tapply(fishing_monitored$nbycatch,fishing_monitored$metiers,length))
}

} #metier else bracket

} # sample iteration loop
if (bymetier==FALSE) {
BPUE_est_mean<-mean(BPUE_est,na.rm=TRUE) #if bymetier is TRUE then 2 element vector
BPUE_est_CV<-sd(BPUE_est,na.rm=TRUE)/mean(BPUE_est,na.rm=TRUE) 
} else {
BPUE_est_mean<-colMeans(BPUE_est,na.rm=TRUE) #if bymetier is TRUE then 2 element vector
BPUE_est_CV<-apply(BPUE_est,2,function(x) sd(x,na.rm=TRUE))/colMeans(BPUE_est,na.rm=TRUE) 
}

return(list(BPUE_est=BPUE_est_mean,CV=BPUE_est_CV))
}




estimate_fishing_effort_metier<-function(fishing=NA,p_report=.9,bymetier=FALSE,metierlevel=TRUE) {
if (metierlevel==FALSE) {
effort<-data.frame(DaS_report=NA,hauls_report=NA,DaS_estimated=NA,hauls_estimated=NA,DaS_real=NA,hauls_real=NA)
#we calculate real effort, reported effort, and estimated effort where we know the fleet size and assume average effort for those unreported
# here reporting bias is at the vessel/boat level at first random then biased by high res metier (eg small vessels not reporting)
if (bymetier==FALSE) {
reporting<-sample(1:nboat,floor(p_report*nboat),replace=FALSE)
fishing.report<-fishing[fishing$boat%in%reporting,]
effort$hauls_real<-dim(fishing)[1]
effort$hauls_report<-dim(fishing.report)[1]
effort$hauls_estimated<-dim(fishing.report)[1]+(nboat-floor(p_report*nboat))*mean(table(fishing.report$boat)) #add 1-p_report nboat fishing mean number of hauls of those reporting
Das<-rowSums(table(fishing$boat,fishing$fishing.day)!=0)
Das_report<-rowSums(table(fishing.report$boat,fishing.report$fishing.day)!=0)

effort$DaS_real<-sum(Das)
effort$DaS_report<-sum(Das_report)
effort$DaS_estimated=sum(Das_report)+(nboat-floor(p_report*nboat))*mean(Das_report)

} else {
#p_report is a vector
nmetier<-length(unique(fishing$metiers))

reporting<-sample(unique(fishing$boat[fishing$metiers==1]),floor(p_report[1]*length(unique(fishing$boat[fishing$metiers==1]))),replace=FALSE) # sample without replacement
for (i in 2:nmetier) {
reporting<-c(reporting,sample(unique(fishing$boat[fishing$metiers==i]),floor(p_report[i]*length(unique(fishing$boat[fishing$metiers==i]))),replace=FALSE)) # sample without replacement

}

#reporting<-sample(1:nboat,floor(p_report*nboat),replace=FALSE)

fishing.report<-fishing[fishing$boat%in%reporting,]
effort$hauls_real<-dim(fishing)[1]
effort$hauls_report<-dim(fishing.report)[1]
effort$hauls_estimated<-dim(fishing.report)[1]+(nboat-floor(p_report*nboat))*mean(table(fishing.report$boat)) #add 1-p_report nboat fishing mean number of hauls of those reporting
Das<-rowSums(table(fishing$boat,fishing$fishing.day)!=0)
Das_report<-rowSums(table(fishing.report$boat,fishing.report$fishing.day)!=0)

effort$DaS_real<-sum(Das)
effort$DaS_report<-sum(Das_report)
effort$DaS_estimated<-sum(Das_report)+(nboat-floor(p_report*nboat))*mean(Das_report)

}
} else { # if metierlevel bracket
effort<-data.frame(metier=unique(fishing$metiers),DaS_report=NA,hauls_report=NA,DaS_estimated=NA,hauls_estimated=NA,DaS_real=NA,hauls_real=NA)
###
if (bymetier==FALSE) {
reporting<-sample(1:nboat,floor(p_report*nboat),replace=FALSE)
fishing.report<-fishing[fishing$boat%in%reporting,]
effort$hauls_real<-as.numeric(table(fishing$metier))
effort$hauls_report<-as.numeric(table(fishing.report$metier))
nvessels_metier<-colSums(table(fishing$boat,fishing$metiers)!=0)-colSums(table(fishing.report$boat,fishing.report$metiers)!=0)
effort$hauls_estimated<-as.numeric(table(fishing.report$metier))+nvessels_metier*colMeans(table(fishing.report$boat,fishing.report$metiers)) #add 1-p_report nboat fishing mean number of hauls of those reporting

#####################!!!
Das<-apply(table(fishing$boat,fishing$fishing.day,fishing$metiers)!=0,c(3,1),sum)
Das_report<-apply(table(fishing.report$boat,fishing.report$fishing.day,fishing.report$metiers)!=0,c(3,1),sum)

effort$DaS_real<-rowSums(Das)
effort$DaS_report<-rowSums(Das_report)
effort$DaS_estimated=rowSums(Das_report)+nvessels_metier*rowMeans(Das_report)

} else {
#p_report is a vector
nmetier<-length(unique(fishing$metiers))

reporting<-sample(unique(fishing$boat[fishing$metiers==1]),floor(p_report[1]*length(unique(fishing$boat[fishing$metiers==1]))),replace=FALSE) # sample without replacement
for (i in 2:nmetier) {
reporting<-c(reporting,sample(unique(fishing$boat[fishing$metiers==i]),floor(p_report[i]*length(unique(fishing$boat[fishing$metiers==i]))),replace=FALSE)) # sample without replacement

}

#reporting<-sample(1:nboat,floor(p_report*nboat),replace=FALSE)

fishing.report<-fishing[fishing$boat%in%reporting,]
effort$hauls_real<-as.numeric(table(fishing$metier))
effort$hauls_report<-as.numeric(table(fishing.report$metier))
nvessels_metier<-colSums(table(fishing$boat,fishing$metiers)!=0)-colSums(table(fishing.report$boat,fishing.report$metiers)!=0)
effort$hauls_estimated<-as.numeric(table(fishing.report$metier))+nvessels_metier*colMeans(table(fishing.report$boat,fishing.report$metiers)) #add 1-p_report nboat fishing mean number of hauls of those reporting

Das<-apply(table(fishing$boat,fishing$fishing.day,fishing$metiers)!=0,c(3,1),sum)
Das_report<-apply(table(fishing.report$boat,fishing.report$fishing.day,fishing.report$metiers)!=0,c(3,1),sum)

effort$DaS_real<-rowSums(Das)
effort$DaS_report<-rowSums(Das_report)
effort$DaS_estimated=rowSums(Das_report)+nvessels_metier*rowMeans(Das_report)

}
###
}

return(effort) #one row if no metier n rows if n metiers
}

##############################################################################3
##############################################################################
#############################################################################

###############################################################
#### simulation code section


p.bycatch<-cbind(c(.1,.01,.001,.001,.001),c(.2,.02,.002,.01,.1))
p.metier<-cbind(c(.5,.8,.2),c(.5,.8,.2))
p_report<-cbind(c(1,.8,.2,.2),c(1,.2,.8,.2))

pmonitor_boatTRUE<-c(.01,.1,.5,.9)
p_monitor_boat_boatTRUE<-c(.01,.1,.5,.9)
pmonitor_boatFALSE<-c(.001,.01,.05,.1,.5,.9)
p_monitor_metier<-cbind(c(1,.8,.2,.2),c(1,.2,.8,.2))


#estimate_fishing_effort_metier<-function(fishing=NA,p_report=.9,bymetier=FALSE,metierlevel=TRUE) {

#monitor_BPUE_metier<-function(pmonitor=0.5,nsample=1000,BPUE_real=0,fishing=NA, p_monitor_boat=.1,boat_samp=TRUE,
#p_haul_obs=1,detect_prob=1,refusal_rate=0, misclassification=0, bymetier=FALSE, p_monitor_metier=.5)

#make_fishing_year_metier<-function(mean.bycatch.event=1,mean.bycatch.large.event=20,p.large.event=0.01,
#							nboat=100,mean.fishing.event.boat.day=2,p.bycatch=c(0.1,.01),p.metier=c(.2,.8),stochastic=TRUE)


#fishing1<-make_fishing_year_metier(p.bycatch=p.bycatch[1,],p.metier=p.metier[1,])

#BPUE_real<-sum(fishing1$nbycatch)/dim(fishing1)[1]

#effort<-estimate_fishing_effort_metier(fishing=fishing1,p_report=p_report[1,],bymetier=TRUE,metierlevel=TRUE)


#monitoring<-monitor_BPUE_metier(pmonitor=pmonitor_boatTRUE[3],nsample=1000,BPUE_real=BPUE_real,fishing=fishing1, 
#p_monitor_boat=p_monitor_boat_boatTRUE[2],boat_samp=TRUE,p_haul_obs=1,detect_prob=1,refusal_rate=0, misclassification=0, bymetier=TRUE, p_monitor_metier=p_monitor_metier[1,])


monitor_estimate<-data.frame(year=NA,p_bycatch_1=NA,p_bycatch_2=NA,p_metier_1=NA,p_metier_2=NA,pmonitor=NA,p_monitor_boat=NA,boat_samp=NA,bymetier=NA,p_monitor_metier=NA,BPUE_real=NA,BPUE_est=NA,BPUE_est_CV=NA)

bymetier<-c("TRUE","FALSE")
boat_samp<-c("TRUE","FALSE")

###
pmonitor_boatFALSE<-c(.9,.5,.1,.05,.01,.001)
###

###
p_monitor_metier<-cbind(c(1,.8,.2,.2),c(1,.2,.8,.2))
###

 #for (y in 1:100) {
#iter<-1
bymetier<-"FALSE"
boat_samp<-"TRUE"
#for (m.s in 1:2) {

#for (b.s in 1:2) {


for (b in 1:dim(p.bycatch)[1]) {

for (m in 1:dim(p.metier)[1]) {
fishing<-make_fishing_year_metier(p.bycatch=p.bycatch[b,],p.metier=p.metier[m,])
BPUE_real<-sum(fishing$nbycatch)/dim(fishing)[1]
print("fishing done")
flush.console()

for (pm.b in 1:length(pmonitor_boatTRUE)) {

for (p_m.b in 1:length(p_monitor_boat_boatTRUE)) {

for (p_m_m in 1:length(p_monitor_metier)) {
if (bymetier==FALSE) {
temp<-data.frame(year=NA,p_bycatch_1=NA,p_bycatch_2=NA,p_metier_1=NA,p_metier_2=NA,pmonitor=NA,p_monitor_boat=NA,boat_samp=NA,bymetier=NA,p_monitor_metier=NA,BPUE_real=NA,BPUE_est=NA,BPUE_est_CV=NA)
} else {
temp<-data.frame(year=rep(NA,2),p_bycatch_1=NA,p_bycatch_2=NA,p_metier_1=NA,p_metier_2=NA,pmonitor=NA,p_monitor_boat=NA,boat_samp=NA,bymetier=NA,p_monitor_metier=NA,BPUE_real=NA,BPUE_est=NA,BPUE_est_CV=NA)

}

temp$p_bycatch_1<-p.bycatch[b,1]
temp$p_bycatch_2<-p.bycatch[b,2]
temp$p_metier_1<-p.metier[m,1]
temp$p_metier_2<-p.metier[m,2]
temp$pmonitor<-pmonitor_boatTRUE[pm.b]
temp$p_monitor_boat<-p_monitor_boat_boatTRUE[p_m.b]
temp$boat_samp<-boat_samp
temp$bymetier<-bymetier
temp$p_monitor_metier<-p_monitor_metier[p_m_m]
temp$BPUE_real<-BPUE_real


monitoring<-monitor_BPUE_metier(pmonitor=pmonitor_boatTRUE[pm.b],nsample=1000,BPUE_real=BPUE_real,fishing=fishing, 
p_monitor_boat=p_monitor_boat_boatTRUE[p_m.b],boat_samp=boat_samp,p_haul_obs=1,detect_prob=1,refusal_rate=0, misclassification=0, bymetier=bymetier, p_monitor_metier=p_monitor_metier[p_m_m,])

temp$BPUE_est<-monitoring$BPUE_est
temp$BPUE_est_CV<-monitoring$CV

monitor_estimate<-rbind(monitor_estimate,temp)

print(b)
print(m)
flush.console()
} #p_m_m
} # p_m.b
} # pm.b
} #m
} #b
#} #b.s
#} #m.s 
 
 #} #y


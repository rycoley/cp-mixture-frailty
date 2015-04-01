##### rycoley@gmail.com
##### Simulation code from Coley and Brown "Latent Class Approach to Compound Poisson Mixture Frailty Model for HIV Prevention Studies""
##### part of dissertation work at University of Washington

##Simulation code for Scenario VI: Estimate single CP frailty model when true data generating mechanism is single CP frailty model, use informative priors centered on true value
##run in parallel with multicore on SGE computing cluster

#this code is analogous to cpe-ds-I.R. Please see annotated version in the simulations folder for this same paper on github for explanation of variable names, functions, etc.


library(multicore)

PNR<-0.65

rho.b<-3
rho.a<-(rho.b*PNR)/(1-PNR) 

lam.s<-0.05*50
lam.r<-50
eta.s<-6 
eta.r<-2 


S<- 500 
A<-10000
B<- 25000 
ke<- 1

source("cpe-source.R")
source("cpe-to-run.R")	###
	
write.csv( "completed", paste("status-",PNR,"-",ER,".csv",sep=""))
	
	
res<-mclapply(1:S,do.one,mc.cores=20)


res.beta<-res.lam<-res.eta<-res.rho<-matrix(nrow=S,ncol=4)
res.HD<-matrix(nrow=S,ncol=3)

for(s in 1:S){
	res.beta[s,]<-res[[s]]$bet.res
	res.lam[s,]<-res[[s]]$lam.res
	res.rho[s,]<-res[[s]]$rho.res
	res.eta[s,]<-res[[s]]$eta.res
	res.HD[s,]<-res[[s]]$HD}


write.csv(res.beta, "beta-cpe-ds-VI.csv") 
write.csv(res.lam, "lam-cpe-ds-VI.csv") 
write.csv(res.eta, "eta-cpe-ds-VI.csv") 
write.csv(res.rho, "rho-cpe-ds-VI.csv") 
write.csv(res.HD, "HD-cpe-ds-VI.csv") 






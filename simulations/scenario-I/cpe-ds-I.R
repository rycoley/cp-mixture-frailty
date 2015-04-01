##### rycoley@gmail.com
##### Simulation code from Coley and Brown "Latent Class Approach to Compound Poisson Mixture Frailty Model for HIV Prevention Studies""
##### part of dissertation work at University of Washington


##Simulation code for Scenario I: Estimate single CP frailty model when true data generating mechanism is latent class CP mixture frailty model, use informative priors centered on true value
##run in parallel with multicore on SGE computing cluster


library(multicore)

##priors
#shape and rate parameters for gamma prior on baseline hazard
lam.s<-0.05*50
lam.r<-50

#shape and rate parameters for gamma prior on eta
eta.s<-6 
eta.r<-2

#beta prior for proportion at risk
#Note: in the paper I put the prior on 1-exp(-rho)=P(Risk). In the code for estimating a single CP frailty distribution, I use the prior for P(No Risk)=exp(-rho). This is merely an artifact in how I originally coded this earlier model. There is no impact on estimates
PR<-0.35
rho.s1<-(3*(1-PR))/PR 
rho.s2<-3


S<- 500  #number of simulations
A<- 10000 #length of burn-in
B<- 25000   #total number of iterations to perform
ke<- 5 #number thinned

library(survival)

source("cpe-ds-source-I.R") #log posterior likelihoods and sampling functions
source("cpe-ds-to-run-I.R")	#do.one() function does data generation and sampling algorithm

res<-mclapply(1:S,do.one,mc.cores=20) #multicore apply function for 20 cores


#matrices in which to store simulation results
res.beta<-res.lam<-res.eta<-res.rho<-matrix(nrow=S, ncol=4)
res.HD<-matrix(nrow=S, ncol=3)


#fill matrices with results from S simulations
for(s in 1:S){
        res.beta[s,] <- res[[s]]$bet.res ###
        res.lam[s,]<-res[[s]]$lam.res
        res.eta[s,]<-res[[s]]$eta.res
        res.rho[s,]<-res[[s]]$rho.res
        res.HD[s,]<-res[[s]]$HD}
		

#write results to csv files
write.csv(res.beta, "beta-cpe-ds-I.csv") 
write.csv(res.lam, "lam-cpe-ds-I.csv") 
write.csv(res.eta, "eta-cpe-ds-I.csv") 
write.csv(res.rho, "rho-cpe-ds-I.csv") 
write.csv(res.HD, "HD-cpe-ds-I.csv") 





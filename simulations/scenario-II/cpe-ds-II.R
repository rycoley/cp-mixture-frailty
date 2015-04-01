##### rycoley@gmail.com
##### Simulation code from Coley and Brown "Latent Class Approach to Compound Poisson Mixture Frailty Model for HIV Prevention Studies""
##### part of dissertation work at University of Washington


##Simulation code for Scenario II: Estimate single CP frailty model when true data generating mechanism is latent class CP mixture frailty model, use noninformative priors
##run in parallel with multicore on SGE computing cluster

#this code is analogous to cpe-ds-I.R. Please see annotated version in the simulations folder for this same paper on github for explanation of variable names, functions, etc.

library(multicore)


rho.s1<-1
rho.s2<-1

S<- 500 
A<- 10000
B<- 25000  
ke<- 5

library(survival)

source("cpe-ds-source-II.R")
source("cpe-ds-to-run-II.R")	###

res<-mclapply(1:S,do.one,mc.cores=20) 

res.beta<-res.lam<-res.eta<-res.rho<-matrix(nrow=S, ncol=4)
res.HD<-matrix(nrow=S, ncol=3)


for(s in 1:S){
        res.beta[s,] <- res[[s]]$bet.res ###
        res.lam[s,]<-res[[s]]$lam.res
        res.eta[s,]<-res[[s]]$eta.res
        res.rho[s,]<-res[[s]]$rho.res
        res.HD[s,]<-res[[s]]$HD}


write.csv(res.beta, "beta-cpe-ds-II.csv") 
write.csv(res.lam, "lam-cpe-ds-II.csv") 
write.csv(res.eta, "eta-cpe-ds-II.csv") 
write.csv(res.rho, "rho-cpe-ds-II.csv") 
write.csv(res.HD, "HD-cpe-ds-II.csv") 

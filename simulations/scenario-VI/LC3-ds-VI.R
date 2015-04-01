##### rycoley@gmail.com
##### Simulation code from Coley and Brown "Latent Class Approach to Compound Poisson Mixture Frailty Model for HIV Prevention Studies""
##### part of dissertation work at University of Washington


##Simulation code for Scenario VI: Estimating 3 latent classes and CP mixture model but shared CP frailty is true data generating mechanism
##run in parallel with multicore on SGE computing cluster

#this code is analogous to LC3-ds-I.R. Please see annotated version in this same folder on github for explanation of variable names, functions, etc.

library(multicore)

mu.gam<-c(log(2),log(1.5), log(3), log(0.5))
sig.gam<-rep(0.5,4)
mu.c<-c(-0.5, 1)
sig.c<-c(1,1)


PR<-c(0.05, 0.25, 0.8)
b<-4
rho.s1<-(b*PR)/(1-PR)
rho.s2<-rep(b,3)

eta.s<-c(4,6,8)
eta.r<-c(2,2,2)

lam.s<-0.05*50
lam.r<-50


S<- 500
A<-10000
B<- 25000 
ke<- 5


source("LC3-ds-source-I.R") #use same source code as scenario I
source("LC3-ds-to-run-VI.R")	


res<-mclapply(1:S,do.one,mc.cores=40)

res.beta<-res.lam<-res.eta1<-res.eta2<-res.eta3<-res.rho1<-res.rho2<-res.rho3<-res.c1<-res.c2<-res.gam1<-res.gam2<-res.gam3<-res.gam4<-res.PC1<-res.PC2<-res.PC3<-matrix(nrow=S, ncol=4)
res.HD<-matrix(nrow=S, ncol=3)


for(s in 1:S){
        res.beta[s,] <- res[[s]]$bet.res ###
        res.lam[s,]<-res[[s]]$lam.res
        res.eta1[s,]<-res[[s]]$eta.res[1,]
        res.eta2[s,]<-res[[s]]$eta.res[2,]
        res.eta3[s,]<-res[[s]]$eta.res[3,]
        res.rho1[s,]<-res[[s]]$rho.res[1,]
        res.rho2[s,]<-res[[s]]$rho.res[2,]
        res.rho3[s,]<-res[[s]]$rho.res[3,]
        res.PC1[s,]<-res[[s]]$PC.res[1,]
        res.PC2[s,]<-res[[s]]$PC.res[2,]
        res.PC3[s,]<-res[[s]]$PC.res[3,]
        res.gam1[s,]<-res[[s]]$gam.res[1,]
        res.gam2[s,]<-res[[s]]$gam.res[2,]
        res.gam3[s,]<-res[[s]]$gam.res[3,]
        res.gam4[s,]<-res[[s]]$gam.res[4,]
        res.c1[s,]<-res[[s]]$c.res[1,]
        res.c2[s,]<-res[[s]]$c.res[2,]
        res.HD[s,]<-res[[s]]$HD}



write.csv(res.beta, "beta-LC3-ds-VI.csv") 
write.csv(res.lam, "lam-LC3-ds-VI.csv") 
write.csv(res.eta1, "eta1-LC3-ds-VI.csv") 
write.csv(res.eta2, "eta2-LC3-ds-VI.csv") 
write.csv(res.eta3, "eta3-LC3-ds-VI.csv") 
write.csv(res.rho1, "rho1-LC3-ds-VI.csv") 
write.csv(res.rho2, "rho2-LC3-ds-VI.csv") 
write.csv(res.rho3, "rho3-LC3-ds-VI.csv") 
write.csv(res.PC1, "PC1-LC3-ds-VI.csv") 
write.csv(res.PC2, "PC2-LC3-ds-VI.csv") 
write.csv(res.PC3, "PC3-LC3-ds-VI.csv") 
write.csv(res.gam1, "gam1-LC3-ds-VI.csv") 
write.csv(res.gam2, "gam2-LC3-ds-VI.csv") 
write.csv(res.gam3, "gam3-LC3-ds-VI.csv") 
write.csv(res.gam4, "gam4-LC3-ds-VI.csv") 
write.csv(res.c1, "c1-LC3-ds-VI.csv") 
write.csv(res.c2, "c2-LC3-ds-VI.csv") 
write.csv(res.HD, "HD-LC3-ds-VI.csv") 


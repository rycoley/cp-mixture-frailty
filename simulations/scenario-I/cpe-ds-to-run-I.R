##### rycoley@gmail.com
##### Simulation code from Coley and Brown "Latent Class Approach to Compound Poisson Mixture Frailty Model for HIV Prevention Studies""
##### part of dissertation work at University of Washington

##Simulation code for Scenario I: Estimate single CP frailty model when true data generating mechanism is latent class CP mixture frailty model, use informative priors centered on true value
#Accompanies cpe-ds-I.R and cpe-ds-source-I.R 
#This file runs a single simulation, including data generation based on random seed and runs gibbs sampler to get parameter estimates

library(compiler)
enableJIT(1)


do.one<-function(seed){#, A, B, ke, rho.s1, rho.s2, lam.s, lam.r, eta.s, eta.r){#

print(c(seed,"assigned")) #for monitoring
set.seed(seed)

### DATA GENERATION
#define true parameter values for data generation and to compare with posterior estimates
lam.true<-lambda<-0.05
bet.true<-beta<-log(0.5)

rho<-c(-log(0.95), -log(0.75), -log(0.2)) #5%, 25%,80% at risk
eta<-c(2,3,4)

gamma<-c(log(2),log(1.5), log(3), log(0.5))
cv<-c(-0.5, 1) #this is referred to as intercept alpha in the paper

#Step 1 in data generation: put participants into latent risk classes
#this matrix gives the probability of membership in each class (column) for each individual(row)
pm<-matrix(nrow=n, ncol=3)
pm[,1]<-expit(cv[1]-as.vector(V%*%gamma))
pm[,2]<-expit(cv[2]-as.vector(V%*%gamma))-expit(cv[1]-as.vector(V%*%gamma))
pm[,3]<-1-expit(cv[2]-as.vector(V%*%gamma))

#this matrix shows which class (column) each individual (row) is assigned to
cm<-t(apply(pm,1,rmultinom,n=1,size=1))

#transform class matrix above into a vector
class.vec<-vector(length=n)
class.vec[cm[,1]==1]<-1
class.vec[cm[,2]==1]<-2
class.vec[cm[,3]==1]<-3

#true value of nu for this simulation run
nu<-(1/n)*(t(as.matrix(cm%*%rho))%*%as.matrix(cm%*%eta))

#Step 2 of data generation: Generate latent frailty terms for each person conditional on her latent class
N<-rep(0,n) #number of risk processes
for(i in 1:n){N[i]<-rpois(1,rho[class.vec[i]])}

Z<-rep(0,n) #sum of risks associated with each risk process
for(i in 1:n){Z[i]<-rgamma(1, shape=N[i]*eta[class.vec[i]], rate=nu)}

#Step 3 of data generation: Generate survival times for each person given her frailty, impose administrative censoring
Y<-rep(NA,n) #survival times
for(i in 1:n){if(!N[i]==0){
		Y[i]<-rexp(1, rate=(lambda*exp(X.mat[i]*beta)*Z[i]))	}
		else {Y[i]<-1000}	}

Tf<-quantile(Y,0.06) #6% event rate
delta<-rep(0,n)
delta[Y<=Tf]<-1		
Y[Y>Tf]<-Tf

#Calculate observed hazard (pdf) using censoring time rounded to months. Save this to compare to simulated hazard
Tf.mo<-c(1:length(months))[months==min(months[Tf<months])]
pdf.true<-get.surv.pdf(cdf=get.surv.mo(Y=get.ymo(Y, Tf.mo=Tf.mo),delta=delta, Tf.mo=Tf.mo), Tf.mo=Tf.mo)

rho.true<-mean(cm%*%rho)
eta.true<-mean(cm%*%eta)

rho<-eta<-lambda<-beta<-N<-Z<-NULL

### GIBB'S SAMPLER
keep<-seq(A,B,ke)
K<-length(keep)

#define matrices in which to save posterior samples and set initial parameter values
rho.mat<-eta.mat<-lambdas<-betas<-vector(length=K)
rhohat<- runif(1,0.25,0.75) 
etahat<- runif(1,2,4) 
lam<- runif(1,0.025,0.075) 
bet<- 0 

#sample latent variables
Nhat<-Zhat<-vector(length=n)
pk0<- woah(rhoi=rhohat, etai=etahat, lexbi=lam, Yi=Tf)
pk1<- woah(rhoi=rhohat, etai=etahat, lexbi=(lam*exp(bet)), Yi=Tf)

for(i in 1:n){ NiZi<-get.NiZi(rhoi=rhohat, etai=etahat, lexbi=lam*exp(X.mat[i]*bet), Yi=Y[i], deltai=delta[i], N0=1)
	Nhat[i]<-NiZi$Nc
	Zhat[i]<-NiZi$Zc}

#vector for saving hellinger distance
HD<-vector(length=K)

#run sampling algorithm
for(j in 2:B){
	
	if(j %in% keep){k<-(j+(ke-1))/ke}
	
	rhohat<-slice.rho(rho0=rhohat, eta=etahat, Nhat=Nhat, Zhat=Zhat, a=rho.s1, b=rho.s2)
	
	etahat<-slice.eta(eta0=etahat, rho=rhohat, Nhat=Nhat, Zhat=Zhat, eta.s=eta.s, eta.r=eta.r)
	
	pk0<- woah(rhoi=rhohat, etai=etahat, lexbi=lam, Yi=Tf)
	pk1<- woah(rhoi=rhohat, etai=etahat, lexbi=(lam*exp(bet)), Yi=Tf)
	
	for(i in 1:n){ NiZi<-get.NiZi(rhoi=rhohat, etai=etahat, lexbi=lam*exp(X.mat[i]*bet), Yi=Y[i], deltai=delta[i], N0=Nhat[i])
	Nhat[i]<-NiZi$Nc
	Zhat[i]<-NiZi$Zc}
	
	lam<-get.lambda(delta=delta, Zs=Zhat, bet=bet, Y=Y, lam.s=lam.s, lam.r=lam.r)
	
	bet<-slice.beta(beta0=bet, lam=lam, Zhat=Zhat, delta=delta, Y=Y)
	
	if(j%in%seq(1,A,ke)){print(j)}
	
	if(j%in%keep){print(j)
		k.num<-c(1:K)[j==seq(A,B,ke)]
		rho.mat[k.num]<-rhohat
		eta.mat[k.num]<-etahat
		lambdas[k.num]<-lam
		betas[k.num]<-bet

#can save latent variable samples as well, but slows down code and not used for sim results
#		N.mat[k.num,]<-Nhat
#		Z.mat[k.num,]<-Zhat

		HD[k.num]<-get.yrep.surv(bet=bet, lam=lam, rhohat=rhohat, etahat=etahat, pdf.true=pdf.true, Tf=Tf, Tf.mo=Tf.mo)
		}
}

#summarize simulation results
bet.res<-c(median(betas), cover(x=betas, true=bet.true))
lam.res<-c(median(lambdas), cover(x=lambdas, true=lam.true))
rho.res<-c(median(rho.mat), cover(x=rho.mat, true=rho.true))
eta.res<-c(median(eta.mat), cover(x=eta.mat, true=eta.true))

HD<-quantile(HD,p=c(0.5,0.025,0.975))

print(c(seed,"completed")) #for monitoring

results<-list(bet.res=bet.res, lam.res=lam.res, rho.res=rho.res, eta.res=eta.res, HD=HD, seed=seed)
return(results)}


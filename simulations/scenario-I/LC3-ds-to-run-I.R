##### rycoley@gmail.com
##### Simulation code from Coley and Brown "Latent Class Approach to Compound Poisson Mixture Frailty Model for HIV Prevention Studies""
##### part of dissertation work at University of Washington

##Simulation code for Scenario I: Correct priors and correctly specified model
#Accompanies LC3-ds-I.R and LC3-ds-source-I.R
#This file runs a single simulation, including data generation based on random seed and runs gibbs sampler to get parameter estimates

library(compiler)
enableJIT(1)
library(survival)


do.one<-function(seed){

print(c(seed,"assigned")) #for monitoring
set.seed(seed)

### DATA GENERATION

#define true parameter values for data generation and to compare with posterior estimates
lam.true<-lambda<-0.05
bet.true<-beta<-log(0.5)
rho.true<-rho<-c(-log(0.95), -log(0.75), -log(0.2)) #5%, 25%,80% at risk
eta.true<-eta<-c(2,3,4)

gam.true<-gamma<-c(log(2),log(1.5), log(3), log(0.5))
c.true<-cv<-c(-0.5, 1) #this is referred to as intercept alpha in the paper

#Step 1 in data generation: put participants into latent risk classes
#this matrix gives the probability of membership in each class (column) for each individual(row)
pm<-matrix(nrow=n, ncol=3)
pm[,1]<-expit(cv[1]-as.vector(V%*%gamma))
pm[,2]<-expit(cv[2]-as.vector(V%*%gamma))-expit(cv[1]-as.vector(V%*%gamma))
pm[,3]<-1-expit(cv[2]-as.vector(V%*%gamma))

#this matrix shows which class (column) each individual (row) is assigned to
cm<-t(apply(pm,1,rmultinom,n=1,size=1))

#proportion in each class for this simulation run
p.class.true<-apply(cm,2,mean)

#transform class matrix above into a vector
class.vec<-vector(length=n)
class.vec[cm[,1]==1]<-1
class.vec[cm[,2]==1]<-2
class.vec[cm[,3]==1]<-3

#true value of nu for this simulation run
nu<-sum(rho.true*eta.true*apply(pm,2,mean))

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


rho<-gamma<-cv<-lambda<-beta<-eta<-nu<-N<-Z<-cm<-class.vec<-NULL


### GIBB'S SAMPLER
keep<-seq(A,B,ke)
K<-length(keep)

#define matrices in which to save posterior samples and set initial parameter values
betas<-lambdas<-vector(length=K)
bet<-0
lam<-runif(1, 0.025, 0.075)

#frailty model parameters
rho.mat<-eta.mat<-p.class<-matrix(nrow=K, ncol=3)
rhohat<-c( -log( runif(1,0.85,1)), -log(runif(1,0.65, 0.85)), -log(runif(1,0.1,0.3)) )
etahat<-c(runif(1,1.5,2.5), runif(1,2.5,3.5), runif(1,3.5,4.5))

#ordinal regression model parameters
gam.mat<-matrix(nrow=K, ncol=D)
gamhat<-runif(D,-0.25, 0.25)

#matrix to save sampled values of intercepts of ordinal regression model, refered to as alpha in paper but with "c" in this code
cs<-matrix(nrow=K, ncol=2)  
cv<-c(runif(1,-1,0), runif(1,0,1)) #cv is the vector of intercepts (alpha_1, alpha_2)

#pm is a matrix that holds mulitnomial probability of classification for each individual (row) and is updated in sampler
pm<-matrix(nrow=n, ncol=3)
pm[,1]<-expit(cv[1]-as.vector(V%*%gamhat))
pm[,2]<-expit(cv[2]-as.vector(V%*%gamhat))-expit(cv[1]-as.vector(V%*%gamhat))
pm[,3]<-1-expit(cv[2]-as.vector(V%*%gamhat))

#cm is a matrix holding the sampled multinomial variable using probability in pm
cm<-t(apply(pm,1,rmultinom,n=1,size=1))

#class.vec transforms cm into a vector
class.vec<-vector(length=n)
class.vec[cm[,1]==1]<-1
class.vec[cm[,2]==1]<-2
class.vec[cm[,3]==1]<-3

#frailty model parameter, must be updated as rho or eta change or classification probabilites change in order to maintain identifiability, E(Z)=1
nu<-sum(rhohat*etahat*apply(pm,2,mean))

#latent variables N,Z
Nhat<-Zhat<-vector(length=n)
for(i in 1:n){
	NiZi<-get.NiZi(rhoi=rhohat[class.vec[i]], etai=etahat[class.vec[i]], nu=nu, lexbi=lam*exp(X.mat[i]*bet), Yi=Y[i], deltai=delta[i], N0=1)
	Nhat[i]<-NiZi$Nc
	Zhat[i]<-NiZi$Zc}

#Hellinger Distance
HD<-vector(length=K)

#run sampling algorithm
for(j in 2:B){

	eta.vec<-as.vector(cm%*%etahat)

	cv<-slice.c(c0=cv, gam=gamhat, class.vec=class.vec, rhoh=rhohat, etah=etahat, eta.vec=eta.vec, Ns=Nhat, Zs=Zhat, mu.c=mu.c, sig.c=sig.c) #vector of intercepts, referred to as alpha in paper

	gamhat<-slice.gam(gam0=gamhat, cv=cv, class.vec=class.vec, rhoh=rhohat, etah=etahat, eta.vec=eta.vec, Ns=Nhat, Zs=Zhat, mu.gam=mu.gam, sig.gam=sig.gam)

	#update classification probabilities and frailty distribution
	pm<-get.pm(cv=cv, gam=gamhat)
	nu<-sum(rhohat*etahat*apply(pm,2,mean))

	#sample classes for each person
	for(i in 1:n){
		class.vec[i]<-MH.Ri(R0=class.vec[i], Vi=V[i,], cv=cv, gam=gamhat, rho=rhohat, Ni=Nhat[i], eta=etahat, nu=nu, Zi=Zhat[i])	}
	for(k in 1:3){cm[,k]<-class.vec==k}

	rhohat<- slice.rho(rho0=rhohat, Ns=Nhat, Zs=Zhat, etah=etahat, class.vec=class.vec, pm=pm, cm=cm, rho.s1=rho.s1, rho.s2=rho.s2)

	etahat<-slice.eta(eta0=etahat, Ns=Nhat, Zs=Zhat, rhoh=rhohat, class.vec=class.vec, pm=pm, cm=cm, eta.s=eta.s, eta.r=eta.r)

	nu<-sum(rhohat*etahat*apply(pm,2,mean))
	
	#sample frailty for each person
for(i in 1:n){
	NiZi<-get.NiZi(rhoi=rhohat[class.vec[i]], etai=etahat[class.vec[i]], nu=nu,  lexbi=lam*exp(X.mat[i]*bet), Yi=Y[i], deltai=delta[i], N0=Nhat[i])
	Nhat[i]<-NiZi$Nc
	Zhat[i]<-NiZi$Zc}

	lam<-rgamma(1,shape=(sum(delta)+lam.s), rate=(sum(Zhat*exp(X.mat*bet)*Y)+lam.r))
		
	bet<-slice.beta(beta0=bet, delta=delta, Nhat=Nhat, Zhat=Zhat, lam=lam, Y=Y)

#		if(j%in%seq(1,A,ke)){print(j)} #for monitoring, if desired

		if(j %in%keep){#print(j)
		k.num<-c(1:K)[j==seq(A,B,ke)]
		cs[k.num,]<-cv
		gam.mat[k.num,]<-gamhat
		rho.mat[k.num,]<-rhohat
		eta.mat[k.num,]<-etahat
		lambdas[k.num]<-lam
		betas[k.num]<-bet

#if you wanted to keep track of class and frailty samples, can slow things down quite a bit
#		R.mat[k.num,]<-class.vec
#		N.mat[k.num,]<-Nhat
#		Z.mat[k.num,]<-round(Zhat,4)

		for(k in 1:3){p.class[k.num,k]<-sum(class.vec==k)/n}
		
		HD[k.num]<-get.yrep.surv(bet=bet, lam=lam, gamhat=gamhat, cv=cv, rhohat=rhohat, etahat=etahat, delta=delta, pdf.true=pdf.true, Tf=Tf, Tf.mo=Tf.mo)}
}

#summarize simulation results
bet.res<-c(median(betas), cover(x=betas, true=bet.true))
lam.res<-c(median(lambdas), cover(x=lambdas, true=lam.true))

rho.res<-eta.res<-PC.res<-matrix(nrow=3,ncol=4)
for(k in 1:3){
	rho.res[k,]<-c(median(rho.mat[,k]), cover(x=rho.mat[,k], true=rho.true[k]))
	eta.res[k,]<-c(median(eta.mat[,k]), cover(x=eta.mat[,k], true=eta.true[k]))
	PC.res[k,]<-c(median(p.class[,k]), cover(x=p.class[,k], true=p.class.true[k]))}


gam.res<-matrix(nrow=D,ncol=4)
for(d in 1:D){gam.res[d,]<-c(median(gam.mat[,d]), cover(x=gam.mat[,d], true=gam.true[d]))}

c.res<-matrix(nrow=2,ncol=4)
for(k in 1:2){c.res[k,]<-c(median(cs[,k]), cover(x=cs[,k], true=c.true[k]))}

HD<-quantile(HD,p=c(0.5,0.025, 0.975))

print(c(seed,"completed"))

results<-list(bet.res=bet.res, lam.res=lam.res, rho.res=rho.res, eta.res=eta.res, PC.res=PC.res, gam.res=gam.res, c.res=c.res, HD=HD, seed=seed) #, mZ=mZ
return(results)}






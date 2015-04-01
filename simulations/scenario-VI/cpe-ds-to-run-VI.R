##### rycoley@gmail.com
##### Simulation code from Coley and Brown "Latent Class Approach to Compound Poisson Mixture Frailty Model for HIV Prevention Studies""
##### part of dissertation work at University of Washington

##Simulation code for Scenario VI: Estimate single CP frailty model when true data generating mechanism is single CP frailty model, use informative priors centered on true value
#Accompanies cpe-ds-VI.R and cpe-ds-source-VI.R
#This file runs a single simulation, including data generation based on random seed and runs gibbs sampler to get parameter estimates

#this code is similar to cpe-ds-to-run-I.R. Please see annotated version in the simulations folder for this same paper on github for explanation of variable names, functions, etc.


library(survival)

do.one<-function(seed){ 

print(c(seed,"assigned")) 
set.seed(seed)

#define parameter values
rho<--log(PNR) 
true.rho<-rho

true.lam<-lambda<-0.05
true.bet<-beta<-log(0.5)

true.eta<-eta<-3
nu<-rho*eta

### GENERATE DATA
N<-rep(0,n) #number of risk processes
for(i in 1:n){N[i]<-rpois(1,rho)}

Z<-rep(0,n) #sum of risks associated with each risk process
for(i in 1:n){Z[i]<-rgamma(1, shape=N[i]*eta, rate=nu)}

Y<-rep(NA,n)
for(i in 1:n){if(!N[i]==0){
		Y[i]<-rexp(1, rate=(lambda*exp(X[i]*beta)*Z[i]))	}
		else {Y[i]<-1000}	}

Tf<-quantile(Y,0.06) #censoring at 6% event rate
delta<-rep(0,n)
delta[Y<=Tf]<-1		
Y[Y>Tf]<-Tf

#we want to compare observed and predicted survival up to Tf.mo
Tf.mo<-c(1:length(months))[months==min(months[Tf<months])]

#hazard function for observed data
pdf.true<-get.surv.pdf(cdf=get.surv.mo(Y=get.ymo(Y, Tf.mo=Tf.mo),delta=delta, Tf.mo=Tf.mo), Tf.mo=Tf.mo)


### GIBB'S SAMPLER

#set-up vectors/matrices for storing sampled values
keep<-seq(A,B,ke)
K<-length(keep)

rho.mat<-eta.mat<-lambdas<-betas<-vector(length=K)
Nhat<-Zhat<-vector(length=n)
HD<-vector(length=K)

rho<-eta<-lambda<-beta<-NULL
N<-Z<-NULL

#get initial values
rhohat<- runif(1,0.5,1) 
etahat<- runif(1,1,3) 
lam<- runif(1,0.025,0.075) 
bet<- 0 

pk0<- woah(rhoi=rhohat, etai=etahat, lexbi=lam, Yi=Tf)
pk1<- woah(rhoi=rhohat, etai=etahat, lexbi=(lam*exp(bet)), Yi=Tf)

	for(i in 1:n){NiZi<-get.NiZi(deltai=delta[i], X=X[i], pk0=pk0, pk1=pk1, rhoi=rhohat, etai=etahat, lexbi=lam*exp(X[i]*bet), Yi=Y[i], N0=1)
		Nhat[i]<-NiZi$Nc
		Zhat[i]<-NiZi$Zc}

for(j in 2:B){
		
	rhohat<-slice.rho(rho0=rhohat, eta=etahat, Nhat=Nhat, Zhat=Zhat, rho.a=rho.a, rho.b=rho.b)

	etahat<-slice.eta(eta0=etahat, rho=rhohat, Nhat=Nhat, Zhat=Zhat, eta.s=eta.s, eta.r=eta.r)
	
	pk0<- woah(rhoi=rhohat, etai=etahat, lexbi=lam, Yi=Tf)
	pk1<- woah(rhoi=rhohat, etai=etahat, lexbi=(lam*exp(bet)), Yi=Tf)
	
	for(i in 1:n){NiZi<-get.NiZi(deltai=delta[i], X=X[i], pk0=pk0, pk1=pk1, rhoi=rhohat, etai=etahat, lexbi=lam*exp(X[i]*bet), Yi=Y[i], N0=Nhat[i])
		Nhat[i]<-NiZi$Nc
		Zhat[i]<-NiZi$Zc}
		
	lam<-get.lambda(delta=delta, Zs=Zhat, bet=bet, Y=Y, lam.s=lam.s, lam.r=lam.r)
	bet<-slice.beta(beta0=bet, lam=lam, Zhat=Zhat, delta=delta, Y=Y)


		if(j %in%keep){#print(j)
		k.num<-c(1:K)[j==seq(A,B,ke)]
		rho.mat[k.num]<-rhohat
		eta.mat[k.num]<-etahat
		lambdas[k.num]<-lam
		betas[k.num]<-bet

		HD[k.num]<-get.yrep.surv(bet=bet, lam=lam, rhohat=rhohat, etahat=etahat,  pdf.true=pdf.true, Tf=Tf, Tf.mo=Tf.mo)		}}

bet.res<-c(median(betas), cover(x=betas, true=true.bet))
lam.res<-c(median(lambdas), cover(x=lambdas, true=true.lam))
rho.res<-c(median(rho.mat), cover(x=rho.mat, true=true.rho))
eta.res<-c(median(eta.mat), cover(x=eta.mat, true=true.eta))

HD<-quantile(HD,p=c(0.5, 0.025, 0.975))


results<-list(bet.res=bet.res, lam.res=lam.res, rho.res=rho.res, eta.res=eta.res, HD=HD, seed=seed) 

write.csv(c(read.csv(paste("status-",PNR,"-",ER,".csv",sep="")),seed), paste("status-",PNR,"-",ER,".csv",sep=""))


return(results)}

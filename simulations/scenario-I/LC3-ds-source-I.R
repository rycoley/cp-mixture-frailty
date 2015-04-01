##### rycoley@gmail.com
##### Simulation code from Coley and Brown "Latent Class Approach to Compound Poisson Mixture Frailty Model for HIV Prevention Studies""
##### part of dissertation work at University of Washington

##Simulation code for Scenario I: Correct priors and correctly specified model
#Accompanies LC3-ds-I.R and LC3-ds-to-run-I.R for correct priors
#Also accompanies LC3-ds-V.R and LC3-ds-to-run-V.R when data is generated with unbalanced classes
#Accompanies LC3-ds-VI.R and LC3-ds-to-run-VI.R when data is generated following single CP frailty distribution


#This file contains the posterior log likelihoods and sampling functions for all model parameters, as well as the sampling functions for latent variables.


#data set-up that doesn't change across simulations
n1<-1500 #number in each arm
n2<-1500 #
n<-n1+n2
X.mat<-c(rep(0,n1), rep(1,n2))	#treatment indicator

#covariates used for classification
age<-rep(c(rep(1,500),rep(2,500), rep(3,500)),2)
sti<-rep(c(rep(1,100),rep(0,400), rep(1,80),rep(0,420), rep(1,60), rep(0,440)),2) 
marr.co<-rep(0,n)
for(i in 1:n){if(age[i]==1){marr.co[i]<-rbinom(1,1,0.25)}
	else{if(age[i]==2){marr.co[i]<-rbinom(1,1,0.5)}
		else{marr.co[i]<-rbinom(1,1,0.75)}}}

V<-as.matrix(cbind(age==1, age==2, sti, marr.co))
D<-dim(V)[2]


###FUNCTIONS, constant across simulations
#defined with fixed values (above) but not simulated data (from LC3-ds-to-run-I.R)

#log-posterior likelihood for log-HR (beta)
#prior(beta)=N(0,10)
lik.beta<-function(bet, h, Ns, Zs, delta, lam, Y){
		XB<-as.vector(X.mat[Ns>0]*bet)
		return( sum(delta[Ns>0]*XB) - sum(Zs[Ns>0] * Y[Ns>0] * lam * exp(XB)) + log(dnorm(bet[h],0,10)) )}

#slice sampler for beta
slice.beta<-function(beta0, delta, lam, Nhat, Zhat, Y){
	#for(h in 1:H){ #would cycle through h for beta with more than one entry
	h<-1
	z<-J<-K<-L<-R<-beta.star<-NULL
	z<- lik.beta(bet=beta0, h=h, Ns=Nhat, Zs=Zhat, delta=delta, lam=lam, Y=Y) - rexp(1,1) #want to make a slice with all betas that have a log-likelihood above this
	w<- 0.05 #width of each step
	m<- 10 #number of total steps out
	betL<-betR<-beta.star<-beta0 
	betL[h]<- beta0[h] - (w*runif(1,0,1)) #starting lower bound of slice
	betR[h]<- betL[h] + w #starting upper bound of slice
	J<- floor(m*runif(1,0,1)) #keep track of number of steps
	K<- (m-1)-J
	while(lik.beta(bet=betL, h=h, Ns=Nhat, Zs=Zhat,  delta=delta, lam=lam, Y=Y) > z & J>0){betL[h]<- betL[h]-w
		J<- J-1} #extend slice to contain all betas below starting beta with log likelihood above z
	while(lik.beta(bet=betR, h=h, Ns=Nhat, Zs=Zhat, delta=delta, lam=lam, Y=Y) > z & K>0){betR[h]<- betR[h]+w
		K<- K-1}
	beta.star[h]<- runif(1,betL[h],betR[h]) #sample new beta from this slice 
	while(lik.beta(bet=beta.star, h=h, Ns=Nhat, Zs=Zhat,  delta=delta, lam=lam, Y=Y) < z){
		if(beta.star[h]<beta0[h]){betL[h]<-beta.star[h]}
		if(beta.star[h]>beta0[h]){betR[h]<-beta.star[h]}
		beta.star[h]<- runif(1,betL[h],betR[h])} #repat sampling to obtain sample of beta that has log-likelihood above z
		beta0<-beta.star#}
	return(beta.star)}

#inverse logit function
expit<-function(x){return(exp(x)/(1+exp(x)))}

#function to populate the matrix containing classification probabilities
get.pm<-function(cv, gam){
	pm<-matrix(nrow=n, ncol=3)
	pm[,1]<-expit(cv[1]-(as.vector(V%*%gam)) )
	pm[,2]<-expit(cv[2]-(as.vector(V%*%gam) ))-expit(cv[1]-(as.vector(V%*%gam)))
	pm[,3]<-1-expit(cv[2]-(as.vector(V%*%gam)))
	return(pm)}

#log likelihood for frailty model parameter nu, which is determined by values of latent variables and other frailty model parameters 
get.lik.nu<-function(nuhat, Ns, Zs, eta.vec){
	return(sum(Ns[Ns>0]*eta.vec[Ns>0]*log(nuhat)) - sum(Zs[Ns>0]*nuhat))}


#log posterior likelihood for intercept terms in regression model (referred to as alpha in paper)
#prior(c)=N(mu.c, sig.c)
lik.c<-function(cv, k, gam, class.vec, rhoh, etah, eta.vec, Ns, Zs, mu.c, sig.c){

	pm<-get.pm(cv=cv, gam=gam)
	nu<-sum(rhoh*etah*apply(pm,2,mean))
	lik.nu<-get.lik.nu(nuhat=nu, Ns=Ns, Zs=Zs, eta.vec=eta.vec)

	if(k==1){return(cv[1]*sum(class.vec==1) - sum(log(1+exp(cv[1]-as.vector(V[class.vec==1,]%*%gam)))) + sum(log(expit(cv[2]-as.vector(V[class.vec==2,]%*%gam)) - expit(cv[1]-as.vector(V[class.vec==2,]%*%gam)))) + lik.nu+ log(dnorm(cv[1], mean=mu.c[1],sd=sig.c[1])) )}
	
	if(k==2){return(sum(log(expit(cv[2]-as.vector(V[class.vec==2,]%*%gam)) - expit(cv[1]-as.vector(V[class.vec==2,]%*%gam)))) + sum(log(1-expit(cv[2]-as.vector(V[class.vec==3,]%*%gam)))) + lik.nu + log(dnorm(cv[2], mean=mu.c[2], sd=sig.c[2])) )}	
		}
	
#slice sampler for intercepts (referred to here as c but as alpha in paper)
slice.c<-function(c0, gam, class.vec, rhoh, etah, eta.vec, Ns, Zs, mu.c, sig.c){
	for(k in 1:2){
	if(k==1){min.c<- -5 
		max.c<- c0[2]-1e-10}
	if(k==2){min.c<- c0[1]+1e-10
		max.c<-5}
	z<-J<-K<-L<-R<-c.star<-NULL
	z<- lik.c(cv=c0, k=k, gam=gam, class.vec=class.vec, rhoh=rhoh, etah=etah, eta.vec=eta.vec, Ns=Ns, Zs=Zs, mu.c=mu.c, sig.c=sig.c) -rexp(1,1) 
	w<- 0.05
	m<- 10
	cL<-cR<-c.star<-c0
	cL[k]<- c0[k] - (w*runif(1,0,1))
	cR[k]<- cL[k] + w
	cL[k]<-max(cL[k], min.c)
	cR[k]<-min(cR[k], max.c)
	J<- floor(m*runif(1,0,1))
	K<- (m-1)-J
	while(cL[k]>min.c & lik.c(cv=cL, k=k, gam=gam, class.vec=class.vec, rhoh=rhoh, etah=etah, eta.vec=eta.vec, Ns=Ns, Zs=Zs, mu.c=mu.c, sig.c=sig.c) > z & J>0){cL[k]<-max(cL[k]-w, min.c)
		J<- J-1}
		cL[k]<-max(cL[k],min.c)
	while(cR[k]<max.c & lik.c(cv=cR, k=k, gam=gam, class.vec=class.vec, rhoh=rhoh, etah=etah, eta.vec=eta.vec, Ns=Ns, Zs=Zs, mu.c=mu.c, sig.c=sig.c) > z & K>0){cR[k]<- min(cR[k]+w, max.c)
		K<- K-1}
		cR[k]<-min(cR[k],max.c)
	c.star[k]<- runif(1,cL[k],cR[k])
	while(lik.c(cv=c.star, k=k, gam=gam, class.vec=class.vec, rhoh=rhoh, etah=etah, eta.vec=eta.vec, Ns=Ns, Zs=Zs, mu.c=mu.c, sig.c=sig.c) < z){
		if(c.star[k] < c0[k]){L[k]<- c.star[k]}
		if(c.star[k] > c0[k]){R[k]<- c.star[k]}
		c.star[k]<- runif(1,cL[k],cR[k])}
		c0<-c.star}
	return(c.star)	}
	
#log posterior likelihood for coefficients gamma in ordinal regression
#prior(gam)=N(mu.gam, sig.gam)
lik.gam<-function(gam, d, cv, class.vec, rhoh, etah, eta.vec, Ns, Zs, mu.gam, sig.gam){
	
	pm<-get.pm(cv=cv, gam=gam)
	nu<-sum(rhoh*etah*apply(pm,2,mean))
	lik.nu<-get.lik.nu(nuhat=nu, Ns=Ns, Zs=Zs, eta.vec=eta.vec)

	return(	sum(-as.vector(V[class.vec==1,]%*%gam)) - sum(log(1+exp(cv[1]-as.vector(V[class.vec==1,]%*%gam)))) + 	sum(log(expit(cv[2]-as.vector(V[class.vec==2,]%*%gam))-expit(cv[1]-as.vector(V[class.vec==2,]%*%gam)))) + sum(log(1-expit(cv[2]-as.vector(V[class.vec==3,]%*%gam)))) + lik.nu + log(dnorm(gam[d],mu.gam, sig.gam)))}


#slice sampler for gamma
slice.gam<-function(gam0, cv, class.vec, rhoh, etah, eta.vec, Ns, Zs, mu.gam, sig.gam){
	for(d in 1:D){
			min.gam<- -2.5
			max.gam<- 2.5
	z<-J<-K<-L<-R<-gam.star<-NULL
	z<- lik.gam(gam=gam0, d=d, cv=cv, class.vec=class.vec, rhoh=rhoh, etah=etah, eta.vec=eta.vec, Ns=Ns, Zs=Zs, mu.gam=mu.gam[d], sig.gam=sig.gam[d]) - rexp(1,1)
	w<- 0.05
	m<- 10
	gamL<-gamR<-gam.star<-gam0
	gamL[d]<- gam0[d] - (w*runif(1,0,1))
	gamR[d]<- gamL[d] + w
	gamL[d]<- max(gamL[d], min.gam)
	gamR[d]<- min(gamR[d], max.gam)
	J<- floor(m*runif(1,0,1))
	K<- (m-1)-J	
	while(gamL[d]>min.gam & lik.gam(gam=gamL, d=d, cv=cv, class.vec=class.vec, rhoh=rhoh, etah=etah, eta.vec=eta.vec, Ns=Ns, Zs=Zs, mu.gam=mu.gam[d], sig.gam=sig.gam[d]) > z & J>0){gamL[d]<- max(gamL[d]-w,min.gam)
		J<- J-1}
		gamL[d]<-max(gamL[d],min.gam)
	while(gamR[d]<max.gam & lik.gam(gam=gamR, d=d,  cv=cv, class.vec=class.vec, rhoh=rhoh, etah=etah, eta.vec=eta.vec, Ns=Ns, Zs=Zs, mu.gam=mu.gam[d], sig.gam=sig.gam[d]) > z & K>0){gamR[d]<- min(gamR[d]+w,max.gam)
		K<- K-1}
		gamR[d]<-min(gamR[d],max.gam)
	gam.star[d]<- runif(1,gamL[d],gamR[d])
	while(lik.gam(gam=gam.star, d=d, cv=cv, class.vec=class.vec, rhoh=rhoh, etah=etah, eta.vec=eta.vec, Ns=Ns, Zs=Zs, mu.gam=mu.gam[d], sig.gam=sig.gam[d]) < z){
		if(gam.star[d]<gam0[d]){gamL[d]<-gam.star[d]}
		if(gam.star[d]>gam0[d]){gamR[d]<-gam.star[d]}
		gam.star[d]<- runif(1,gamL[d],gamR[d])}
		gam0<-gam.star}
	return(gam.star)}


#log posterior likelihood for mean number of exposure processes within each class
#prior(PR)=beta(A,B)
lik.rho<-function(rho, k, Ns, Zs, etah, class.vec, pm, cm, rho.s1, rho.s2){
	nu<-sum(rho*etah*apply(pm,2,mean))
	return(sum(Ns[class.vec==k]*log(rho[k])) - rho[k]*sum(class.vec==k) + sum(Ns[Ns>0]*as.vector(cm[Ns>0,]%*%etah) *log(nu) ) - sum(Zs[Ns>0]*nu) + (rho.s1-1)*log(1-exp(-rho[k])) - rho.s2*rho[k]) } #using prior put on probability at risk

#slice sampler for rho
slice.rho<-function(rho0, Ns, Zs, etah, class.vec, pm, cm, rho.s1, rho.s2){
	for(k in 1:3){
		if(k==1){min.rho<-1e-20
			max.rho<-rho0[2]}
		if(k==2){min.rho<- rho0[1]
				max.rho<-rho0[3]}
		if(k==3){min.rho<-rho0[2]
			max.rho<-5}		
	z<-J<-K<-L<-R<-rho.star<-NULL
	z<- lik.rho(rho=rho0, k=k, Ns=Ns, Zs=Zs, etah=etah, class.vec=class.vec, pm=pm, cm=cm, rho.s1=rho.s1[k], rho.s2=rho.s2[k]) - rexp(1,1)
	w<- 0.05
	m<- 10
	rhoL<-rhoR<-rho.star<-rho0
	rhoL[k]<- rho0[k] - (w*runif(1,0,1))
	rhoR[k]<- rhoL[k] + w
	rhoL[k]<- max(rhoL[k], min.rho)
	rhoR[k]<- min(rhoR[k], max.rho)
	J<- floor(m*runif(1,0,1))
	K<- (m-1)-J	
	while(rhoL[k]>min.rho & lik.rho(rho=rhoL, k=k, Ns=Ns, Zs=Zs, etah=etah, class.vec=class.vec, pm=pm, cm=cm, rho.s1=rho.s1[k], rho.s2=rho.s2[k])> z & J>0){rhoL[k]<- max(rhoL[k]-w,min.rho)
		J<- J-1}
		rhoL[k]<-max(rhoL[k],min.rho)
	while(rhoR[k]<max.rho & lik.rho(rho=rhoR, k=k, Ns=Ns, Zs=Zs, etah=etah, class.vec=class.vec, pm=pm, cm=cm, rho.s1=rho.s1[k], rho.s2=rho.s2[k]) > z & K>0){rhoR[k]<- min(rhoR[k]+w,max.rho)
		K<- K-1}
		rhoR[k]<-min(rhoR[k],max.rho)
	rho.star[k]<- runif(1,rhoL[k],rhoR[k])
	while(lik.rho(rho=rho.star, k=k, Ns=Ns, Zs=Zs, etah=etah, class.vec=class.vec, pm=pm, cm=cm, rho.s1=rho.s1[k], rho.s2=rho.s2[k]) < z){
		if(rho.star[k]<rho0[k]){rhoL[k]<-rho.star[k]}
		if(rho.star[k]>rho0[k]){rhoR[k]<-rho.star[k]}
		rho.star[k]<- runif(1,rhoL[k],rhoR[k])}
		rho0<-rho.star}
	return(rho.star)}

#log posterior likelihood for eta, the shape paramter of the gamma random variable for amount of risk associated with each exposure process within latent class
#prior(eta)=Gamma(O_{eta}, T_{eta}) = Gamma(eta.s, eta.r)
lik.eta<-function(eta, k, rhoh, Ns, Zs, pm, cm, class.vec, eta.s, eta.r){
	nu<-sum(rhoh*eta*apply(pm,2,mean))
	return(sum(Ns[Ns>0 & class.vec==k]*eta[k]*log(Zs[Ns>0 & class.vec==k]+1e-100)) -sum(lgamma(Ns[Ns>0 & class.vec==k]*eta[k])) + sum(Ns[Ns>0]*as.vector(cm[Ns>0,]%*%eta)*log(nu)) -sum(Zs[Zs>0]*nu) + log(dgamma(eta[k], shape=eta.s, rate=eta.r)))}

#slice sampler for eta
slice.eta<-function(eta0, Ns, Zs, rhoh, class.vec, pm, cm, eta.s, eta.r){
	for(k in 1:3){
		if(k==1){min.eta<-1e-20
			max.eta<-eta0[2]}
		if(k==2){min.eta<- eta0[1]
				max.eta<-eta0[3]}
		if(k==3){min.eta<-eta0[2]
			max.eta<-100}		
	z<-J<-K<-L<-R<-eta.star<-NULL
	z<- lik.eta(eta=eta0, k=k,rhoh=rhoh,  Ns=Ns, Zs=Zs, pm=pm, cm=cm, class.vec=class.vec, eta.s=eta.s[k], eta.r=eta.r[k]) - rexp(1,1)
	w<- 0.05
	m<- 10
	etaL<-etaR<-eta.star<-eta0
	etaL[k]<- eta0[k] - (w*runif(1,0,1))
	etaR[k]<- etaL[k] + w
	etaL[k]<- max(etaL[k], min.eta)
	etaR[k]<- min(etaR[k], max.eta)
	J<- floor(m*runif(1,0,1))
	K<- (m-1)-J	
	while(etaL[k]>min.eta & lik.eta(eta=etaL, k=k,rhoh=rhoh,  Ns=Ns, Zs=Zs, pm=pm, cm=cm, class.vec=class.vec, eta.s=eta.s[k], eta.r=eta.r[k])> z & J>0){etaL[k]<- max(etaL[k]-w,min.eta)
		J<- J-1}
		etaL[k]<-max(etaL[k],min.eta)
	while(etaR[k]<max.eta & lik.eta(eta=etaR, k=k,rhoh=rhoh,  Ns=Ns, Zs=Zs, pm=pm, cm=cm, class.vec=class.vec, eta.s=eta.s[k], eta.r=eta.r[k]) > z & K>0){etaR[k]<- min(etaR[k]+w,max.eta)
		K<- K-1}
		etaR[k]<-min(etaR[k],max.eta)
	eta.star[k]<- runif(1,etaL[k],etaR[k])
	while(lik.eta(eta=eta.star,k=k,rhoh=rhoh,  Ns=Ns, Zs=Zs, pm=pm, cm=cm, class.vec=class.vec, eta.s=eta.s[k], eta.r=eta.r[k]) < z){
		if(eta.star[k]<eta0[k]){etaL[k]<-eta.star[k]}
		if(eta.star[k]>eta0[k]){etaR[k]<-eta.star[k]}
		eta.star[k]<- runif(1,etaL[k],etaR[k])}
		eta0<-eta.star}
	return(eta.star)}

#log-transformed density for the number of exposure processes for participant i, marginalized over Z[i]
lik.Ni<-function(Ni, etai, deltai, rhoi,  nu, lexbi, Yi){
	return(Ni*(log(rhoi) + etai*log(nu)) + deltai*log(Ni*etai) - (Ni*etai + deltai)*log(nu + lexbi*Yi) -lgamma(Ni+1) )}	

#proposal probabilities for metropolis-hastings algorithm
pN<-function(given){
	if (given==1){return(1/2)}
	if (given>1){return(1/3)}}

#metropolis-hastings algorithm to sample N[i], the number of exposure processes for participant i
MH.Ni<-function(N0, etai, deltai, rhoi, nu, lexbi, Yi){
	U<-runif(1,0,1)
	if(N0>1){
		if(U < 1/3 ){Nc<- N0-1 }
		if(U > 1/3 & U < 2/3){Nc<- N0
			return(Nc)}
		if(U > 2/3){Nc<- N0+1}}
	if(N0==1){
		if(U < 1/2){Nc<- N0
			return(Nc)}
		if(U > 1/2){Nc<- N0+1}}
	
	accept<- rbinom(1,1, min( (exp(lik.Ni(Ni=Nc, etai=etai, deltai=deltai, rhoi=rhoi, nu=nu, lexbi=lexbi, Yi=Yi )-lik.Ni(Ni=N0, etai=etai, deltai=deltai, rhoi=rhoi, nu=nu, lexbi=lexbi, Yi=Yi ) )) * (pN(Nc)/pN(N0)) ,1) )
	if(accept==1){return(Nc)}
	if(accept==0){return(N0) }} 

#mean parameter for Poisson dist on N[i] when no event observed, delta[i]=0
woah<-function(rhoi, nu, etai, lexbi, Yi){ 
	return(rhoi*(nu/(nu+(lexbi*Yi)))^(etai))} 

#function to sample Ni and Zi
get.NiZi<-function(rhoi,etai,nu,lexbi,Yi,deltai,N0){
	if(deltai==0){ #only those without event observed
		Nc<-rpois(1, woah(rhoi=rhoi, nu=nu, etai=etai, lexbi=lexbi, Yi=Yi))	
		if(Nc==0){Zc<-0}}
		if(deltai==1){
		mh.res<- MH.Ni(N0=N0, etai=etai, deltai=deltai, rhoi=rhoi, nu=nu, lexbi=lexbi, Yi=Yi)		
		Nc<- mh.res[1]		}
	if(Nc>0){Zc<- rgamma(1, shape=(Nc*etai + deltai), rate=(nu + lexbi*Yi))
		if(Zc<1e-100){Zc<-1e-100}}
	return(list(Nc=Nc,Zc=Zc))}	 

#go through similar procedure (as with N above) to sample class for individual i (Ri)
#log-transformed density of latent class 
lik.Ri<-function(Ri, Vi, cv, gam, rho, Ni, eta, nu, Zi){
	if(Ri==1){ log.pRi<- cv[1]- Vi%*%gam - log(1+exp(cv[1]-Vi%*%gam)) }
	if(Ri==2){log.pRi<- log(expit(cv[2]-Vi%*%gam) - expit(cv[1]-Vi%*%gam)) }
	if(Ri==3){log.pRi<- log(1-expit(cv[2]-Vi%*%gam)) }
	if(Ni==0){return(log.pRi - rho[Ri])}
	if(Ni>0 & Zi>1e-100){return(log.pRi + Ni*log(rho[Ri]) - rho[Ri] + (Ni*eta[Ri]*(log(nu)+ log(Zi))) -lgamma(Ni*eta[Ri]) )}
	if(Ni>0 & Zi<1e-100){return(log.pRi + Ni*log(rho[Ri]) - rho[Ri] + (Ni*eta[Ri]*(log(nu)+ log(Zi+1e-100))) -lgamma(Ni*eta[Ri]) )}}

#proposal density
pR<-function(given){
	if(given==2){return(1/3)}
	else{return(1/2)}}

#Metropolis-Hastings algorithm for sampling class
MH.Ri<-function(R0, Vi, cv, gam, rho, Ni, eta, nu, Zi){ 
	U<-runif(1,0,1)
	
	if(R0==1){if(U<1/2){return(R0)}
		else{Rc<-2}}
	if(R0==2){if(U<1/3){Rc<-1}
		else{if(U<2/3){return(R0)}
			else{Rc<-3}}}
	if(R0==3){if(U<1/2){Rc<-2}
		else{return(R0)}}
		
	accept<-rbinom(1,1,min( exp(lik.Ri(Ri=Rc, Vi=Vi, cv=cv, gam=gam, rho=rho, Ni=Ni, eta=eta, nu=nu, Zi=Zi) - lik.Ri(Ri=R0, Vi=Vi, cv=cv, gam=gam, rho=rho, Ni=Ni, eta=eta, nu=nu, Zi=Zi)) *(pR(Rc)/pR(R0)) , 1 )) 
	if(accept==1){return(Rc)}
	if(accept==0){return(R0)}}



months<-seq(1/12,10,1/12)
months<-round(months,4)

#round event times to months within a year
get.ymo<-function(Y, Tf.mo){
	Y.mo<-rep((1/12),n)
	for(j in 2:Tf.mo){
		Y.mo[Y>months[(j-1)] & Y<=months[j]]<-months[j]}
	Y.mo[Y>months[Tf.mo]]<-months[Tf.mo]
	return(round(Y.mo,4))}

#kaplan meier survival estimate for event times (in months)
get.surv.mo<-function(Y.mo,delta, Tf.mo){
	so<-Surv(time=Y.mo, event=delta)
	km<-survfit(formula=so~1)
	if(length(km$surv)==Tf.mo){return(km$surv)}
	else{km.surv<-vector(length=Tf.mo)
		km$time<-round(km$time,4)
		for(j in 1:Tf.mo){if(sum(km$time==months[j])>0){km.surv[j]<-km$surv[km$time==months[j]]}
			else{if(j==1){km.surv[j]<-1}
				else{km.surv[j]<-km.surv[(j-1)]}}}
				return(km.surv)}}

#calculate density for discrete time from survival function
get.surv.pdf<-function(cdf, Tf.mo){
	pdf<-rep(0,(Tf.mo+1))
	pdf[1]<-1-cdf[1]
	for(t in 2:Tf.mo){
		pdf[t]<-cdf[(t-1)]-cdf[t]} 
	pdf[(Tf.mo+1)]<-cdf[Tf.mo]
	return(pdf)}

#calculate hellinger distance
get.HD<-function(p,q){ 
		hd<-(1/sqrt(2)) * sqrt( sum( (sqrt(p)-sqrt(q))^2 ) )
		return(hd)}
		
#simulate survival times given current parameter samples, censoring
#calculated hellinger distance between observerd survival and simulated  	
get.yrep.surv<-function(bet, lam, gamhat, cv, rhohat, etahat, delta, pdf.true, Tf, Tf.mo){
	R.rep<-N.rep<-Z.rep<-rep(0,n)
	
	pm<-matrix(nrow=n, ncol=3)
	pm[,1]<-expit(cv[1]-as.vector(V%*%gamhat))
	pm[,2]<-expit(cv[2]-as.vector(V%*%gamhat))-expit(cv[1]-as.vector(V%*%gamhat))
	pm[,3]<-1-expit(cv[2]-as.vector(V%*%gamhat))

	cm<-t(apply(pm,1,rmultinom,n=1,size=1))
	
	p.class.true<-apply(cm,2,mean)

	R.rep<-vector(length=n)
	R.rep[cm[,1]==1]<-1
	R.rep[cm[,2]==1]<-2
	R.rep[cm[,3]==1]<-3

	rho.vec<-as.vector(cm%*%rhohat)
	eta.vec<-as.vector(cm%*%etahat)
	nu<-sum(rhohat*etahat*apply(pm,2,mean))
	
	N.rep<-rpois(n,rho.vec)
	for(i in 1:n){if(N.rep[i]>0){Z.rep[i]<-rgamma(1,shape=N.rep[i]*eta.vec[i], rate=nu)}}
	
	Y.rep<-rep(500,n)
	for(i in 1:n){if(Z.rep[i]>0){Y.rep[i]<-rexp(1, rate=lam*Z.rep[i]*exp(X.mat[i]*bet))}}

	delta.rep<-rep(0,n)
	delta.rep[Y.rep<=Tf]<-1		
	Y.rep[Y.rep>Tf]<-Tf
	
	HD<-get.HD(p=pdf.true, q=get.surv.pdf(cdf=get.surv.mo(Y.mo=get.ymo(Y.rep, Tf.mo=Tf.mo),delta=delta.rep, Tf.mo=Tf.mo), Tf.mo=Tf.mo))
			
	return(HD=HD)}

#does the posterior sampling distribution cover the true value?
cover<-function(x,true){
	return(c(quantile(x,p=0.025), quantile(x,p=0.975), quantile(x,p=0.025)<true & quantile(x,p=0.975)>true))}


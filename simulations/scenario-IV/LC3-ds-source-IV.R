##### rycoley@gmail.com
##### Simulation code from Coley and Brown "Latent Class Approach to Compound Poisson Mixture Frailty Model for HIV Prevention Studies""
##### part of dissertation work at University of Washington

##Simulation code for Scenario IV: Estimate model with two classes (true data generating mechanism is still 3 classes)
#Accompanies LC3-ds-IV.R and LC3-ds-to-run-IV.R
#This file contains the posterior log likelihoods and sampling functions for all model parameters, as well as the sampling functions for latent variables.

#this code is analogous to LC3-ds-source-I.R. Please see annotated version in the simulations folder for this same paper on github for explanation of variable names, functions, etc.

n1<-1500 
n2<-1500 
n<-n1+n2
X.mat<-c(rep(0,n1), rep(1,n2))	

age<-rep(c(rep(1,500),rep(2,500), rep(3,500)),2)

sti<-rep(c(rep(1,100),rep(0,400), rep(1,80),rep(0,420), rep(1,60), rep(0,440)),2) 

marr.co<-rep(0,n)
for(i in 1:n){if(age[i]==1){marr.co[i]<-rbinom(1,1,0.25)}
	else{if(age[i]==2){marr.co[i]<-rbinom(1,1,0.5)}
		else{marr.co[i]<-rbinom(1,1,0.75)}}}

V<-as.matrix(cbind(age==1, age==2, sti, marr.co))
D<-dim(V)[2]


###FUNCTIONS

lik.beta<-function(bet, h, Ns, Zs, delta, lam, Y){
		XB<-as.vector(X.mat[Ns>0]*bet)
		return( sum(delta[Ns>0]*XB) - sum(Zs[Ns>0] * Y[Ns>0] * lam * exp(XB)) + log(dnorm(bet[h],0,10)) )}


slice.beta<-function(beta0, delta, lam, Nhat, Zhat, Y){
	#for(h in 1:H){
	h<-1
	z<-J<-K<-L<-R<-beta.star<-NULL
	z<- lik.beta(bet=beta0, h=h, Ns=Nhat, Zs=Zhat, delta=delta, lam=lam, Y=Y) - rexp(1,1)
	w<- 0.05
	m<- 10
	betL<-betR<-beta.star<-beta0
	betL[h]<- beta0[h] - (w*runif(1,0,1))
	betR[h]<- betL[h] + w
	J<- floor(m*runif(1,0,1))
	K<- (m-1)-J
	while(lik.beta(bet=betL, h=h, Ns=Nhat, Zs=Zhat,  delta=delta, lam=lam, Y=Y) > z & J>0){betL[h]<- betL[h]-w
		J<- J-1}
	while(lik.beta(bet=betR, h=h, Ns=Nhat, Zs=Zhat, delta=delta, lam=lam, Y=Y) > z & K>0){betR[h]<- betR[h]+w
		K<- K-1}
	beta.star[h]<- runif(1,betL[h],betR[h])
	while(lik.beta(bet=beta.star, h=h, Ns=Nhat, Zs=Zhat,  delta=delta, lam=lam, Y=Y) < z){
		if(beta.star[h]<beta0[h]){betL[h]<-beta.star[h]}
		if(beta.star[h]>beta0[h]){betR[h]<-beta.star[h]}
		beta.star[h]<- runif(1,betL[h],betR[h])}
		beta0<-beta.star#}
	return(beta.star)}

expit<-function(x){return(exp(x)/(1+exp(x)))}


get.pm<-function(cv, gam){
	pm<-matrix(nrow=n, ncol=2)
	pm[,1]<-expit(cv[1]-(as.vector(V%*%gam)) )
	pm[,2]<-1-pm[,1]
	return(pm)}

get.lik.nu<-function(nuhat, Ns, Zs, eta.vec){
	return(sum(Ns[Ns>0]*eta.vec[Ns>0]*log(nuhat)) - sum(Zs[Ns>0]*nuhat))}


lik.c<-function(cv, gam, class.vec, rhoh, etah, eta.vec, Ns, Zs, mu.c, sig.c){

	pm<-get.pm(cv=cv, gam=gam)
	nu<-sum(rhoh*etah*apply(pm,2,mean))
	lik.nu<-get.lik.nu(nuhat=nu, Ns=Ns, Zs=Zs, eta.vec=eta.vec)

	return(cv*sum(class.vec==1) - sum(log(1+exp(cv-as.vector(V[class.vec==1,]%*%gam)))) + sum(log(1 - expit(cv-as.vector(V[class.vec==2,]%*%gam)))) + lik.nu + log(dnorm(cv, mean=mu.c,sd=sig.c)) )}
	
slice.c<-function(c0, gam, class.vec, rhoh, etah, eta.vec, Ns, Zs, mu.c, sig.c){
	z<-J<-K<-L<-R<-c.star<-NULL
	z<- lik.c(cv=c0, gam=gam, class.vec=class.vec, rhoh=rhoh, etah=etah, eta.vec=eta.vec, Ns=Ns, Zs=Zs, mu.c=mu.c, sig.c=sig.c) -rexp(1,1) 
	w<- 0.05
	m<- 10
	cL<-cR<-c.star<-c0
	cL<- c0 - (w*runif(1,0,1))
	cR<- cL + w
	J<- floor(m*runif(1,0,1))
	K<- (m-1)-J
	while( lik.c(cv=cL, gam=gam, class.vec=class.vec, rhoh=rhoh, etah=etah, eta.vec=eta.vec, Ns=Ns, Zs=Zs, mu.c=mu.c, sig.c=sig.c) > z & J>0){cL<-cL-w
		J<- J-1}
	while(lik.c(cv=cR, gam=gam, class.vec=class.vec, rhoh=rhoh, etah=etah, eta.vec=eta.vec, Ns=Ns, Zs=Zs, mu.c=mu.c, sig.c=sig.c) > z & K>0){cR<- cR+w
		K<- K-1}
	c.star<- runif(1,cL,cR)
	while(lik.c(cv=c.star, gam=gam, class.vec=class.vec, rhoh=rhoh, etah=etah, eta.vec=eta.vec, Ns=Ns, Zs=Zs, mu.c=mu.c, sig.c=sig.c) < z){
		if(c.star < c0){L<- c.star}
		if(c.star > c0){R<- c.star}
		c.star<- runif(1,cL,cR)}
	return(c.star)	}
	

lik.gam<-function(gam, d, cv, class.vec, rhoh, etah, eta.vec, Ns, Zs, mu.gam, sig.gam){
	
	pm<-get.pm(cv=cv, gam=gam)
	nu<-sum(rhoh*etah*apply(pm,2,mean))
	lik.nu<-get.lik.nu(nuhat=nu, Ns=Ns, Zs=Zs, eta.vec=eta.vec)

	return(	sum(-as.vector(V[class.vec==1,]%*%gam)) - sum(log(1+exp(cv-as.vector(V[class.vec==1,]%*%gam)))) + sum(log(1-expit(cv-as.vector(V[class.vec==2,]%*%gam)))) + lik.nu + log(dnorm(gam[d],mu.gam, sig.gam)))}


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


lik.rho<-function(rho, k, Ns, Zs, etah, class.vec, pm, cm, rho.s1, rho.s2){
	nu<-sum(rho*etah*apply(pm,2,mean))
	return(sum(Ns[class.vec==k]*log(rho[k])) - rho[k]*sum(class.vec==k) + sum(Ns[Ns>0]*as.vector(cm[Ns>0,]%*%etah) *log(nu) ) - sum(Zs[Ns>0]*nu) + (rho.s1-1)*log(1-exp(-rho[k])) - rho.s2*rho[k]) } #this is for probability at risk


slice.rho<-function(rho0, Ns, Zs, etah, class.vec, pm, cm, rho.s1, rho.s2){
	for(k in 1:2){
		if(k==1){min.rho<-1e-20
			max.rho<-rho0[2]}
		if(k==2){min.rho<-rho0[1]
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

lik.eta<-function(eta, k, rhoh, Ns, Zs, pm, cm, class.vec, eta.s, eta.r){
	nu<-sum(rhoh*eta*apply(pm,2,mean))
	return(sum(Ns[Ns>0 & class.vec==k]*eta[k]*log(Zs[Ns>0 & class.vec==k]+1e-100)) -sum(lgamma(Ns[Ns>0 & class.vec==k]*eta[k])) + sum(Ns[Ns>0]*as.vector(cm[Ns>0,]%*%eta)*log(nu)) -sum(Zs[Zs>0]*nu) + log(dgamma(eta[k], shape=eta.s, rate=eta.r)))}

slice.eta<-function(eta0, Ns, Zs, rhoh, class.vec, pm, cm, eta.s, eta.r){
	for(k in 1:2){
		if(k==1){min.eta<-1e-20
			max.eta<-eta0[2]}
		if(k==2){min.eta<-eta0[1]
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

lik.Ni<-function(Ni, etai, deltai, rhoi,  nu, lexbi, Yi){
	return(Ni*(log(rhoi) + etai*log(nu)) + deltai*log(Ni*etai) - (Ni*etai + deltai)*log(nu + lexbi*Yi) -lgamma(Ni+1) )}	

pN<-function(given){
	if (given==1){return(1/2)}
	if (given>1){return(1/3)}}

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


woah<-function(rhoi, nu, etai, lexbi, Yi){ #poisson dist for delta=0
	return(rhoi*(nu/(nu+(lexbi*Yi)))^(etai))} 

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

lik.Ri<-function(Ri, Vi, cv, gam, rho, Ni, eta, nu, Zi){
	if(Ri==1){ log.pRi<- cv- Vi%*%gam - log(1+exp(cv-Vi%*%gam)) }
	if(Ri==2){log.pRi<- 1 - expit(cv-Vi%*%gam) }
	if(Ni==0){return(log.pRi - rho[Ri])}
	if(Ni>0 & Zi>1e-100){return(log.pRi + Ni*log(rho[Ri]) - rho[Ri] + (Ni*eta[Ri]*(log(nu)+ log(Zi))) -lgamma(Ni*eta[Ri]) )}
	if(Ni>0 & Zi<1e-100){return(log.pRi + Ni*log(rho[Ri]) - rho[Ri] + (Ni*eta[Ri]*(log(nu)+ log(Zi+1e-100))) -lgamma(Ni*eta[Ri]) )}}


pR<-function(given){
	if(given==2){return(1/3)}
	else{return(1/2)}}

MH.Ri<-function(R0, Vi, cv, gam, rho, Ni, eta, nu, Zi){
	U<-runif(1,0,1)
	
	if(R0==1){if(U<1/2){return(R0)}
		else{Rc<-2}}
	if(R0==2){if(U<1/2){Rc<-1}
		else{return(R0)}}
		
	accept<-rbinom(1,1,min( exp(lik.Ri(Ri=Rc, Vi=Vi, cv=cv, gam=gam, rho=rho, Ni=Ni, eta=eta, nu=nu, Zi=Zi) - lik.Ri(Ri=R0, Vi=Vi, cv=cv, gam=gam, rho=rho, Ni=Ni, eta=eta, nu=nu, Zi=Zi)) *(pR(Rc)/pR(R0)) , 1 )) 
	if(accept==1){return(Rc)}
	if(accept==0){return(R0)}}

months<-seq(1/12,10,1/12)
months<-round(months,4)

get.ymo<-function(Y, Tf.mo){
	Y.mo<-rep((1/12),n)
	for(j in 2:Tf.mo){
		Y.mo[Y>months[(j-1)] & Y<=months[j]]<-months[j]}
	Y.mo[Y>months[Tf.mo]]<-months[Tf.mo]
	return(round(Y.mo,4))}

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

get.surv.pdf<-function(cdf, Tf.mo){
	pdf<-rep(0,(Tf.mo+1))
	pdf[1]<-1-cdf[1]
	for(t in 2:Tf.mo){
		pdf[t]<-cdf[(t-1)]-cdf[t]} 
	pdf[(Tf.mo+1)]<-cdf[Tf.mo]
	return(pdf)}

get.HD<-function(p,q){
		hd<-(1/sqrt(2)) * sqrt( sum( (sqrt(p)-sqrt(q))^2 ) )
		return(hd)}


get.yrep.surv<-function(bet, lam, gamhat, cv, rhohat, etahat, pdf.true, Tf, Tf.mo){
	R.rep<-N.rep<-Z.rep<-rep(0,n)
	
	pm<-matrix(nrow=n, ncol=2)
	pm[,1]<-expit(cv-as.vector(V%*%gamhat))
	pm[,2]<-1-expit(cv-as.vector(V%*%gamhat))

	cm<-t(apply(pm,1,rmultinom,n=1,size=1))
	
	p.class.true<-apply(cm,2,mean)

	R.rep<-vector(length=n)
	R.rep[cm[,1]==1]<-1
	R.rep[cm[,2]==1]<-2

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


cover<-function(x,true){
	return(c(quantile(x,p=0.025), quantile(x,p=0.975), quantile(x,p=0.025)<true & quantile(x,p=0.975)>true))}


##### rycoley@gmail.com
##### Simulation code from Coley and Brown "Latent Class Approach to Compound Poisson Mixture Frailty Model for HIV Prevention Studies""
##### part of dissertation work at University of Washington

##Simulation code for Scenario VI: Estimate single CP frailty model when true data generating mechanism is single CP frailty model, use informative priors centered on true value
#Accompanies cpe-ds-VI.R and cpe-ds-to-run-VI.R 
#This file contains the posterior log likelihoods and sampling functions for all model parameters, as well as the sampling functions for latent variables.

#this code is analogous to cpe-ds-source-I.R. Please see annotated version in the simulations folder for this same paper on github for explanation of variable names, functions, etc.


n1<-1500 
n2<-1500 
n<-n1+n2
X<-c(rep(0,n1), rep(1,n2)) 	

####### FUNCTIONS, also constant across simualations
#defined with fixed values (above) but not simulated data 

lik.beta<-function(bet, lam, Zs, delta, Y){
	XB<-X[Zs>0]*bet
	return(sum(delta[Zs>0]*XB) - sum(Zs[Zs>0]*exp(XB)*lam*Y[Zs>0]) + log(dnorm(bet,0,10)) )}

slice.beta<-function(beta0, lam, Zhat, delta, Y){ 
	z<-J<-K<-L<-R<-beta.star<-NULL
	z<- lik.beta(bet=beta0, lam=lam, Zs=Zhat, delta=delta, Y=Y) - rexp(1,1)
	w<- 0.05
	m<- 10
	L<- beta0 - (w*runif(1,0,1))
	R<- L + w
	J<- floor(m*runif(1,0,1))
	K<- (m-1)-J
	while(lik.beta(bet=L,  lam=lam, Zs=Zhat, delta=delta, Y=Y) > z & J>0){L<- L-w
		J<- J-1}
	while(lik.beta(bet=R,  lam=lam, Zs=Zhat, delta=delta, Y=Y) > z & K>0){R<- R+w
		K<- K-1}
	beta.star<- runif(1,L,R)
	while(lik.beta(bet=beta.star, lam=lam, Zs=Zhat, delta=delta, Y=Y) < z){
		if(beta.star<beta0){L<-beta.star}
		if(beta.star>beta0){R<-beta.star}
		beta.star<- runif(1,L,R)}
	return(beta.star)}

get.lambda<-function(delta, Zs, bet, Y, lam.s, lam.r){
	return(rgamma(1,shape=sum(delta)+lam.s, rate=sum(Zs*exp(X*bet)*Y)+lam.r ))}


ldrho<-function(rhoh,rho.a,rho.b){ 
	return(-rho.a*rhoh + (rho.b-1)*log(1-exp(-rhoh)))} ##this is for probability not at risk! different parameterization than given in paper


lik.rho<-function(rhoh, eta, Ns, Zs, rho.a, rho.b){
	nu<-rhoh*eta
	return(-(n*rhoh)+ sum(Ns[Ns>0]*log(rhoh)) + sum(Ns[Ns>0]*eta*log(nu))   -sum(Zs[Zs>0]*nu) + ldrho(rhoh=rhoh, rho.a=rho.a, rho.b=rho.b) )}


slice.rho<-function(rho0, eta, Nhat, Zhat, rho.a, rho.b){#delta, nuhat, lexbT, 
	z<-J<-K<-L<-R<-rho.star<-NULL
	z<- lik.rho(rhoh=rho0, eta=eta, Ns=Nhat, Zs=Zhat, rho.a=rho.a, rho.b=rho.b) -rexp(1,1) 	
	w<- 0.1
	m<- 10
	L<- rho0 - (w*runif(1,0,1))
	R<- L + w
	L<- max(L,1e-20)
	J<- floor(m*runif(1,0,1))
	K<- (m-1)-J
	while(L>1e-20 & lik.rho(rhoh=L,  eta=eta, Ns=Nhat, Zs=Zhat, rho.a=rho.a, rho.b=rho.b) > z & J>0){
		L<- max(L-w,1e-20)
		J<- J-1}
		L<- max(L, 1e-20)
	while(lik.rho(rhoh=R,  eta=eta, Ns=Nhat, Zs=Zhat, rho.a=rho.a, rho.b=rho.b) > z & R<5 & K>0){R<- R+w
		K<- K-1}
		R<-min(R,5)
	rho.star<- runif(1,L,R)
	while(lik.rho(rhoh=rho.star,  eta=eta, Ns=Nhat, Zs=Zhat, rho.a=rho.a, rho.b=rho.b) < z){
		if(rho.star < rho0){L<- rho.star}
		if(rho.star > rho0){R<- rho.star}
		rho.star<- runif(1,L,R)}
	return(rho.star)	}

lik.eta<-function(etah, rho, Ns, Zs, eta.s, eta.r){
	nu<-rho*etah
	return(sum(Ns[Ns>0]*etah*log(nu)) - sum(lgamma(Ns[Ns>0]*etah)) + sum(Ns[Zs>0]*etah*log(Zs[Zs>0])) -sum(Zs[Zs>0]*nu) + log(dgamma(etah, shape=eta.s, rate=eta.r)) )}


slice.eta<-function(eta0, rho, Nhat, Zhat, eta.s, eta.r){
	z<-J<-K<-L<-R<-eta.star<-NULL
	z<- lik.eta(etah=eta0, rho=rho, Ns=Nhat, Zs=Zhat, eta.s=eta.s, eta.r=eta.r) -rexp(1,1) 	
	w<- 0.2
	m<- 10
	L<- eta0 - (w*runif(1,0,1))
	R<- L + w
	L<-max(L,1e-20)
	J<- floor(m*runif(1,0,1))
	K<- (m-1)-J
	while(L>1e-20 & lik.eta(etah=L, rho=rho, Ns=Nhat, Zs=Zhat, eta.s=eta.s, eta.r=eta.r) > z  & J>0){L<- max(L-w,1e-20)
		J<- J-1}
		L<- max(L, 1e-20)
	while(lik.eta(etah=R, rho=rho, Ns=Nhat, Zs=Zhat, eta.s=eta.s, eta.r=eta.r) > z & K>0){R<- R+w
		K<- K-1}
	eta.star<- runif(1,L,R)
	while(lik.eta(etah=eta.star, rho=rho, Ns=Nhat, Zs=Zhat, eta.s=eta.s, eta.r=eta.r) < z){
		if(eta.star < eta0){L<- eta.star}
		if(eta.star > eta0){R<- eta.star}
		eta.star<- runif(1,L,R)}
	return(eta.star)	}

lik.Ni<-function(Ni,rho,eta,lexbi,Yi){ #for deltai=1
	nu<-rho*eta
	return(Ni*log(rho) -lgamma(Ni+1) + log(Ni*eta) + Ni*eta*log(nu) - ((Ni*eta)+1)*log(nu + (lexbi*Yi))  ) 	}


pN<-function(given){
	if (given==1){return(1/2)}
	if (given>1){return(1/3)}}

MH.Ni<-function(N0, rhoi, etai, lexbi, Yi){
	U<-runif(1,0,1)
	if(N0>1){
		if(U < 1/3 ){Nc<- N0-1 }
		if(U > 1/3 & U < 2/3){Nc<- N0
			return(c(Nc,1))}
		if(U > 2/3){Nc<- N0+1}}
	if(N0==1){
		if(U < 1/2){Nc<- N0
			return(c(Nc,1))}
		if(U > 1/2){Nc<- N0+1}}
	
	accept<- rbinom(1,1, min( (exp(lik.Ni(Ni=Nc,rho=rhoi, eta=etai, lexbi=lexbi, Yi=Yi)-lik.Ni(Ni=N0, rho=rhoi, eta=etai, lexbi=lexbi, Yi=Yi) )) * (pN(Nc)/pN(N0)) ,1) )
	if(accept==1){return(c(Nc,1))}
	if(accept==0){return(c(N0,0)) }}


woah<-function(rhoi, etai, lexbi, Yi){ #poisson dist for delta=0
	nu<-rhoi*etai
	return(rhoi*(nu/(nu+lexbi*Yi))^etai)} 

get.NiZi<-function(deltai, Xi, pk0, pk1, rhoi, etai, lexbi, Yi, N0){
	if(deltai==0){ #only those without event observed
		if(Xi==0){Nc<-rpois(1,pk0)}	
		if(Xi==1){Nc<-rpois(1,pk1)}
		if(Nc==0){Zc<-0}}

	if(deltai==1){
		mh.res<- MH.Ni(N0=N0, rhoi=rhoi, etai=etai, lexbi=lexbi,Yi=Yi)		
		Nc<- mh.res[1]		}

	if(Nc>0){Zc<- rgamma(1, shape=(Nc*etai + deltai), rate=(rhoi*etai + lexbi*Yi))
		if(Zc<1e-100){Zc<-1e-100}}
	return(list(Nc=Nc,Zc=Zc))}	 		


months<-seq(1/12,10,1/12)
months<-round(months,4)

get.ymo<-function(Y, Tf.mo){
	Y.mo<-rep((1/12),n)
	for(j in 2:Tf.mo){
		Y.mo[Y>months[(j-1)] & Y<=months[j]]<-months[j]}
	Y.mo[Y>months[Tf.mo]]<-months[Tf.mo]
	return(Y.mo)}

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
		

get.yrep.surv<-function(bet, lam, rhohat, etahat, ER, pdf.true, Tf, Tf.mo){
	N.rep<-Z.rep<-rep(0,n)
	N.rep<-rpois(n,rhohat)
	nu<-rhohat*etahat
	for(i in 1:n){if(N.rep[i]>0){Z.rep[i]<-rgamma(1,shape=N.rep[i]*etahat, rate=nu)}}
	
	Y.rep<-rep(500,n)
	for(i in 1:n){if(Z.rep[i]>1e-20){Y.rep[i]<-rexp(1, rate=lam*Z.rep[i]*exp(X[i]*bet))}
		if(is.na(Y.rep[i])){Y.rep[i]<-5000}}

	delta.rep<-rep(0,n)
	delta.rep[Y.rep<=Tf]<-1		
	Y.rep[Y.rep>Tf]<-Tf
	
	HD<-get.HD(p=pdf.true, q=get.surv.pdf(cdf=get.surv.mo(Y.mo=get.ymo(Y.rep, Tf.mo=Tf.mo),delta=delta.rep, Tf.mo=Tf.mo), Tf.mo=Tf.mo))
	
	return(HD)}

cover<-function(x,true){
	return(c(quantile(x,p=0.025), quantile(x,p=0.975), quantile(x,p=0.025)<true & quantile(x,p=0.975)>true))}
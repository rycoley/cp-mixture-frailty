##### rycoley@gmail.com
##### Simulation code from Coley and Brown "Latent Class Approach to Compound Poisson Mixture Frailty Model for HIV Prevention Studies""
##### part of dissertation work at University of Washington

##Simulation code for Scenario II: Estimate single CP frailty model when true data generating mechanism is latent class CP mixture frailty model, use noninformative priors
#Accompanies cpe-ds-II.R and cpe-ds-source-II.R 
#This file runs a single simulation, including data generation based on random seed and runs gibbs sampler to get parameter estimates

#this code is analogous to cpe-ds-to-run-I.R. Please see annotated version in the simulations folder for this same paper on github for explanation of variable names, functions, etc.

library(compiler)
enableJIT(1)


do.one<-function(seed){
	
print(c(seed,"assigned"))
set.seed(seed)

#DATA GENERATION
lam.true<-lambda<-0.05
bet.true<-beta<-log(0.5)

rho<-c(-log(0.95), -log(0.75), -log(0.2)) #5%, 25%,80% at risk
eta<-c(2,3,4)

gamma<-c(log(2),log(1.5), log(3), log(0.5))
cv<-c(-0.5, 1)

pm<-matrix(nrow=n, ncol=3)
pm[,1]<-expit(cv[1]-as.vector(V%*%gamma))
pm[,2]<-expit(cv[2]-as.vector(V%*%gamma))-expit(cv[1]-as.vector(V%*%gamma))
pm[,3]<-1-expit(cv[2]-as.vector(V%*%gamma))

cm<-t(apply(pm,1,rmultinom,n=1,size=1))

class.vec<-vector(length=n)
class.vec[cm[,1]==1]<-1
class.vec[cm[,2]==1]<-2
class.vec[cm[,3]==1]<-3

nu<-(1/n)*(t(as.matrix(cm%*%rho))%*%as.matrix(cm%*%eta))

N<-rep(0,n) #number of risk processes
for(i in 1:n){N[i]<-rpois(1,rho[class.vec[i]])}

Z<-rep(0,n) #sum of risks associated with each risk process
for(i in 1:n){Z[i]<-rgamma(1, shape=N[i]*eta[class.vec[i]], rate=nu)}

Y<-rep(NA,n)
for(i in 1:n){if(!N[i]==0){
		Y[i]<-rexp(1, rate=(lambda*exp(X.mat[i]*beta)*Z[i]))	}
		else {Y[i]<-1000}	}

Tf<-quantile(Y,0.06) 
delta<-rep(0,n)
delta[Y<=Tf]<-1		
Y[Y>Tf]<-Tf

Tf.mo<-c(1:length(months))[months==min(months[Tf<months])]
pdf.true<-get.surv.pdf(cdf=get.surv.mo(Y=get.ymo(Y, Tf.mo=Tf.mo),delta=delta, Tf.mo=Tf.mo), Tf.mo=Tf.mo)

rho.true<-mean(cm%*%rho)
eta.true<-mean(cm%*%eta)

rho<-eta<-lambda<-beta<-N<-Z<-NULL

### GIBB'S SAMPLER
keep<-seq(A,B,ke)
K<-length(keep)

rho.mat<-eta.mat<-lambdas<-betas<-vector(length=K)
rhohat<- runif(1,0.25,0.75) 
etahat<- runif(1,2,4) 
lam<- runif(1,0.025,0.075) 
bet<- 0 

Nhat<-Zhat<-vector(length=n)

pk0<- woah(rhoi=rhohat, etai=etahat, lexbi=lam, Yi=Tf)
pk1<- woah(rhoi=rhohat, etai=etahat, lexbi=(lam*exp(bet)), Yi=Tf)

for(i in 1:n){ NiZi<-get.NiZi(rhoi=rhohat, etai=etahat, lexbi=lam*exp(X.mat[i]*bet), Yi=Y[i], deltai=delta[i], N0=1)
	Nhat[i]<-NiZi$Nc
	Zhat[i]<-NiZi$Zc}

HD<-vector(length=K)

for(j in 2:B){
	
	if(j %in% keep){k<-(j+(ke-1))/ke}
	
	rhohat<-slice.rho(rho0=rhohat, eta=etahat, Nhat=Nhat, Zhat=Zhat, a=rho.s1, b=rho.s2)
	
	etahat<-slice.eta(eta0=etahat, rho=rhohat, Nhat=Nhat, Zhat=Zhat)
	
	pk0<- woah(rhoi=rhohat, etai=etahat, lexbi=lam, Yi=Tf)
	pk1<- woah(rhoi=rhohat, etai=etahat, lexbi=(lam*exp(bet)), Yi=Tf)
	
	for(i in 1:n){ NiZi<-get.NiZi(rhoi=rhohat, etai=etahat, lexbi=lam*exp(X.mat[i]*bet), Yi=Y[i], deltai=delta[i], N0=Nhat[i])
	Nhat[i]<-NiZi$Nc
	Zhat[i]<-NiZi$Zc}
	
	lam<-get.lambda(delta=delta, Zs=Zhat, bet=bet, Y=Y)
	
	bet<-slice.beta(beta0=bet, lam=lam, Zhat=Zhat, delta=delta, Y=Y)
	
	if(j%in%seq(1,A,ke)){print(j)}
	
	if(j%in%keep){print(j)
		k.num<-c(1:K)[j==seq(A,B,ke)]
		rho.mat[k.num]<-rhohat
		eta.mat[k.num]<-etahat
		lambdas[k.num]<-lam
		betas[k.num]<-bet

		HD[k.num]<-get.yrep.surv(bet=bet, lam=lam, rhohat=rhohat, etahat=etahat, pdf.true=pdf.true, Tf=Tf, Tf.mo=Tf.mo)

		}
}


bet.res<-c(median(betas), cover(x=betas, true=bet.true))
lam.res<-c(median(lambdas), cover(x=lambdas, true=lam.true))
rho.res<-c(median(rho.mat), cover(x=rho.mat, true=rho.true))
eta.res<-c(median(eta.mat), cover(x=eta.mat, true=eta.true))

HD<-quantile(HD,p=c(0.5,0.025,0.975))

print(c(seed,"completed"))

results<-list(bet.res=bet.res, lam.res=lam.res, rho.res=rho.res, eta.res=eta.res, HD=HD, seed=seed) 
return(results)}

CRT.counts.fit<-function(y=1,t=1,d=1,L=1,coords=1,cluster=1,pair=1,A=1,surround.type="depth",radius=1,prior.cl="weak",prior.sp="weak",npost=1000,buffer_dist=0.1,standard=FALSE){

# returns independent samples from the marginal posterior of all parameters of the extended model
  
# y              an integer vector with the observed counts
# t              a numeric (0/1) vector  with the intervention indicator (1=intervention, 0=control) 
# coords         the matrix with 2 columns for the spatial point locations  
# L              an integer vector of denominators (exposures) (L is at least one) 
# cluster        the numeric vector or factor indicating cluster membership
# pair           the numeric vector or factor indicating pair membership  
# surround type  the surroundedness measure used. Two types only: 'depth' for Tukey's halfspace depth and 'disc'  for  the number of points within a certain radius
# radius         the radius for the disc surroundedness measure
# npost          the number of independent posterior samples 
# prior.cl       the informativeness of the prior for the variance of the cluster random effects. 'weak', 'medium' or 'strong' as per specification of marginal std dev of random effect
# prior.sp       the informativeness of the prior for the variance of the spatial random effects.
  # Optional    
# A is the adjacency (neighbours) matrix for the spatial points 
# d an integer vector with values of a surroundedness measure computed externally
# buffer_dist    the buffer distance for voronoi calculations. See ?st_buffer  

qx<-function(x){
  quantile(x,probs=c(0.5,0.025,0.975))
}
  
n<-length(y)
id<-1:n

# surroundedness measure calculation    
d<-surround(coords,t,type=surround.type,radius=radius)  
d.contra<-(1-t)*d
d.ipsi<-t*d 


tnb<-tess_nb(coords,nb.type="queen",buffer_dist=buffer_dist)
tess<-tnb$tess
A<-tnb$A


# prior setting
msd.weak=0.3 # small since is the standard deviation of the (marginal) random effect (otherwise explodes)
msd.medium=0.05 # target prior marginal standard deviation for the random effect (after integrating the hyperparameter)
msd.strong=0.01 # target marginal standard deviation for the random effect (after integrating the hyperparameter)
# the pc prior paper suggests 0.3 for no good reason

# function only valid when alpha=0.01 (for other values see maple sheet pc_priors_calculations.mw)
u.weak = msd.weak/0.31 
u.medium = msd.medium/0.31 
u.strong = msd.strong/0.31

# a smaller value of msd will imply a stronger penalty for deviating from null model (no random effects)
alpha = 0.01

pc.weak       <-list(prec  = list(prior = "pc.prec", param = c(u.weak, alpha)))
pc.medium     <-list(prec  = list(prior = "pc.prec", param = c(u.medium, alpha)))
pc.strong     <-list(prec  = list(prior = "pc.prec", param = c(u.strong, alpha)))

if (prior.cl=="weak"){
  prior.cl<-pc.weak
}else{
  if (prior.cl=="medium"){
    prior.cl<-pc.medium
  }else{
    prior.cl<-pc.strong
  }
}

if (prior.sp=="weak"){
  prior.sp<-pc.weak
}else{
  if (prior.sp=="medium"){
    prior.sp<-pc.medium
  }else{
    prior.sp<-pc.strong
  }
}


if (standard){
  
  X<-cbind(1,t)
  HQ.orth.std<-orthCV2(X,A)
  
  # model fitting
  dat<-tibble(RESPONSE=y,INTERV=factor(t),CLUSTER=cluster,PAIR=pair,L=L,id=id)
  
  eq.std     <-RESPONSE~INTERV+f(CLUSTER,model="iid",hyper=prior.cl) +f(PAIR,model="iid",hyper=prior.cl)
  eq.std.spat<-RESPONSE~INTERV+f(CLUSTER,model="iid",hyper=prior.cl) +f(PAIR,model="iid",hyper=prior.cl)+ f(id,model="z",Z=HQ.orth.std,hyper=prior.sp)
  
  inla.std      <-inla(eq.std,     data=dat,family="poisson",E=L,control.compute=list(config=TRUE,dic=T,cpo=T))
  inla.std.spat <-inla(eq.std.spat,data=dat,family="poisson",E=L,control.compute=list(config=TRUE,dic=T,cpo=T))
  
  # posterior sampling
  ips.std<-inla.posterior.sample(npost,inla.std)
  
  alpha.samp.std        <-rep(0,npost)
  beta.int.samp.std     <-rep(0,npost)
  
  for (j in 1:npost){
    alpha.samp.std[j]   <-ips.std[[j]]$latent["(Intercept):1",]
    beta.int.samp.std[j]<-ips.std[[j]]$latent["INTERV1:1",]
  }
  
  ips.std.spat<-inla.posterior.sample(npost,inla.std.spat)
  
  alpha.samp.spat        <-rep(0,npost)
  beta.int.samp.spat     <-rep(0,npost)
  
  for (j in 1:npost){
    alpha.samp.spat[j]   <-ips.std.spat[[j]]$latent["(Intercept):1",]
    beta.int.samp.spat[j]<-ips.std.spat[[j]]$latent["INTERV1:1",]
  }
  
  
  coeffs<-data.frame(median=1,q_0.025=1,q_0.975=1)
  coeffs[1,1:3]<-qx(beta.int.samp.std)
  coeffs[2,1:3]<-qx(alpha.samp.std)
  coeffs[3,1:3]<-qx( beta.int.samp.spat)
  coeffs[4,1:3]<-qx(alpha.samp.spat)
  
  rownames(coeffs)<-c("Total effect","Expected count isolated","Total effect (spatial)","Expected count isolated (spatial)")
  
  print(coeffs)
  
  res<-list(coeffs=coeffs,tau.samp=beta.int.samp.std,tau.spat.samp=beta.int.samp.spat,alpha.samp=alpha.samp.std,alpha.spat.samp=alpha.samp.spat,A=A,d=d,tess=tess)

  
}else{
   
  X.ie<-cbind(1,t,d.ipsi,d.contra)
  avg_surr<-avg.surround(t,d)   

# orthogonality algorithmn
HQ.orth<-orthCV2(X.ie,A)


# model fitting
dat<-tibble(RESPONSE=y,INTERV=factor(t),contra=d.contra,ipsi=d.ipsi,CLUSTER=cluster,PAIR=pair,id=id,L=L)
eq <-RESPONSE~INTERV+ipsi+contra+f(CLUSTER,model="iid",hyper=prior.cl) +f(PAIR,model="iid",hyper=prior.cl)+ f(id,model="z",Z=HQ.orth,hyper=prior.sp)
inla.extended <-inla(eq,data=dat,family="poisson",E=L,control.compute=list(config=TRUE,dic=T,cpo=T))

# posterior sampling
ips<-inla.posterior.sample(npost,inla.extended)

alpha.samp        <-rep(0,npost)
beta.int.samp    <-rep(0,npost)
beta.contra.samp <-rep(0,npost)
beta.ipsi.samp   <-rep(0,npost)
tau.samp       <-rep(0,npost)
kappa.samp       <-rep(0,npost)
xi_0.samp        <-rep(0,npost)
xi_1.samp        <-rep(0,npost)


for (j in 1:npost){
  alpha.samp[j]<-ips[[j]]$latent["(Intercept):1",]
  beta.int.samp[j]<-ips[[j]]$latent["INTERV1:1",]
  beta.contra.samp[j]<-ips[[j]]$latent["contra:1",]
  beta.ipsi.samp[j]<-ips[[j]]$latent["ipsi:1",]
  tau.samp[j]<-tot.interv.extended(c(beta.int.samp[j], beta.ipsi.samp[j], beta.contra.samp[j]),t=t,d=d,L=L)
  kappa.samp[j]<-tau.samp[j]-beta.int.samp[j]
  xi_1.samp[j]<- beta.ipsi.samp[j]*avg_surr$ipsi
  xi_0.samp[j]<- beta.contra.samp[j]*avg_surr$contra
}



coeffs<-data.frame(median=1,q_0.025=1,q_0.975=1)
coeffs[1,1:3]<-qx(tau.samp)
coeffs[2,1:3]<-qx(beta.ipsi.samp)
coeffs[3,1:3]<-qx(beta.contra.samp)
coeffs[4,1:3]<-qx(xi_0.samp)
coeffs[5,1:3]<-qx(xi_1.samp)
coeffs[6,1:3]<-qx(beta.int.samp)
coeffs[7,1:3]<-qx(kappa.samp)
coeffs[8,1:3]<-qx(alpha.samp)

rownames(coeffs)<-c("Total effect","Ipsilateral (pair)","Contralateral (pair)","Ipsilateral (average)","Contralateral (average)","Total effect (isolated)","Residual indirect","Expected count isolated")

print(coeffs)

res<-list(coeffs=coeffs,beta.samp=beta.int.samp,eta.samp=beta.ipsi.samp,gamma.samp=beta.contra.samp,tau.samp=tau.samp,kappa.samp=kappa.samp,xi0.samp=xi_0.samp,xi1.samp=xi_1.samp,alpha.samp=alpha.samp,A=A,d=d,tess=tess)

}

  return(res)
}
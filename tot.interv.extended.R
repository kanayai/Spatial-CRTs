tot.interv.extended<-function(theta,t,d,L){
  # for   given parameter values
  # returns the total intervention effect for the extended model
  # and for a  given surroundedness metric d
  
  beta<-theta[1]
  eta<-theta[2]
  gama<-theta[3]
  
  ind.int<-which(t==1)
  ind.cont<-which(t==0)
  d.contra<-d*(1-t)
  d.ipsi<-d*t
  # note alpha cancels out as well as any constant variance
  mu<-exp(beta*t+gama*d.contra+eta*d.ipsi)
  res<-log((sum(mu[ind.int])/sum(L[ind.int]))/(sum(mu[ind.cont])/sum(L[ind.cont])))
  res
}

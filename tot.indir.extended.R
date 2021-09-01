tot.indir.extended<-function(eta=0,gamma=0,t,d){
  # for  a given pairwise indirect effects 'eta' and 'gamma' 
  # returns a list with the two total indirect effects 
  # ipsi and contra
  # for the extended model and for a  given surroundedness metric d
  # 
  
  
  ind.int<-which(t==1)
  ind.cont<-which(t==0)
  
  diff.int<-matrix(NA,n.interv,n.interv)
  diff.cont<-matrix(NA,n.interv,n.interv)
  
  for (i in 2:n.interv){
    for (j  in 1:(i-1)){
      diff.int[i,j] <-abs(d[ind.int[i]]-d[ind.int[j]])
      
      diff.cont[i,j]<-abs(d[ind.cont[i]]-d[ind.cont[j]])
      
    }
  }
  
  ipsi<-eta*mean(diff.int,na.rm=TRUE)
  
  contra<-gamma*mean(diff.cont,na.rm=TRUE)
  
  return(list(ipsi=ipsi,contra=contra))
  
}

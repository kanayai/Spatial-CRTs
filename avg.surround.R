avg.surround<-function(t,d){
  # returns a list with the 2 average absolute difference of pairs 
  #  for a  given surroundedness metric d
  #  one for intervention  and one for control
  
  ind.int<-which(t==1)
  ind.cont<-which(t==0)
  n.interv<-length(ind.int)
  
  diff.int<-matrix(NA,n.interv,n.interv)
  diff.cont<-matrix(NA,n.interv,n.interv)
  
  for (i in 2:n.interv){
    for (j  in 1:(i-1)){
      diff.int[i,j] <-abs(d[ind.int[i]]-d[ind.int[j]])
      
      diff.cont[i,j]<-abs(d[ind.cont[i]]-d[ind.cont[j]])
      
    }
  }
  
  ipsi<-mean(diff.int,na.rm=TRUE)
  contra<-mean(diff.cont,na.rm=TRUE)
  
  return(list(ipsi=ipsi,contra=contra))
  
}


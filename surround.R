surround<-function(X,t,type="depth",radius=1){
  
  n<-length(t)  
  ind.int<-which(t==1)
  n.interv<-length(ind.int)
  d<-rep(0,n)
  
  if (type=="depth"){
    
    for (i in 1:n){
      d[i]<-depth(c(X[i,1],X[i,2]), x=X[setdiff(ind.int,i),], method = "Tukey", approx = FALSE)
      
      if (t[i]==0){
        d[i]<-d[i]*n.interv
      }else{
        d[i]<-d[i]*(n.interv-1) 
      }
    }
    
  }else{
    
    dd<-as.matrix(dist(X))
    
    for (i in 1:n){
      
      if (t[i]==0){
        d[i]<-length(which(dd[i,ind.int]<=radius))
      }
      if (t[i]==1){
        d[i]<-(length(which(dd[i,ind.int]<=radius))-1)
      }
    }
    
  }
  
  return(d)
  
}

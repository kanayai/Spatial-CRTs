orthCV2<-function(ref.mat,A=1,reduce=TRUE,lambda=0){
  # Returns an orthogonalised matrix with respect to a reference matrix
  # A is the adjacency matrix to build the matrix to be orthogonalised
  # ref.mat is the matrix we want to orthogonalised with respect to
  # convergence control parameters
  n.it<-250
  n.ite.moran<-20
  tol<-1e-5
  ite.moran<-1
  ############################ construction of the orthogonal projection matrix
  TS<-ref.mat
  r<-as.numeric(rankMatrix(TS)) # rank of column space
  H<-TS%*%solve(t(TS)%*%TS)%*%t(TS) 
  n<-dim(H)[1]
  # orthogonal complement
  Proj<-diag(n)-H
  r<-n-r # rank of orthogonal complement
  #### option to reduce dimension or not
  if (reduce){
    Moran.op <-Proj%*%A%*%Proj
  }else{
    Moran.op <-Proj%*%Proj%*%Proj
  }
 
  # Construction of the matrix H_{ICAR}
    
  UU<-matrix(0,n,n)
  diag(UU)<-A%*%matrix(1,n,1)
  Q<-UU-A # precision matrix for the ICAR process
  eigQ<-eigen(Q)
  EQ<-eigQ$vectors
  lQ<-eigQ$values
  rr<-qr(Q)$rank # this only keeps non zero eigenvalues
  EQ<-EQ[,1:rr]
  lQ<-lQ[1:rr]
  
  H<-EQ%*%diag(1/sqrt(lQ)) 
  

  ## Construction of the basis matrix for the orthogonal complement
  L<-eigen(Proj)$vectors[,1:r]
 #L<-eigen(Proj)$vectors[,1:r]
  
  ll<-1
  
  while((ll>0)&(ite.moran<n.ite.moran)){
  ite<-1
  nor<-1
  while ((nor>tol)&(ite<n.it)){
     
     # orthogonalisation

      aux<-tcrossprod(L)%*%tcrossprod(H)%*%tcrossprod(L)
      EE<-eigen(aux)
      rr<-length(which(EE$values>1e-10))
      H<-EE$vectors[,1:rr]%*%diag(sqrt(EE$values[1:rr])) 
  
       # the CV property
  
      DD<-diag(1/sqrt(diag(tcrossprod(H))))
      H<-crossprod(DD,H)
  
      nor<-max(abs(crossprod(H,ref.mat)))
      print(nor)
      ite<-ite+1 
  }
  ## Moran reduction
  
  moran.vec<-diag(t(H)%*%Moran.op%*%H)/diag(t(H)%*%Proj%*%H)
  
  ll<-length(which(moran.vec<=lambda))
  # If positive then convert to Moran property
  H<-H[,which(moran.vec>lambda)]
  print(ll)
  ite.moran<-ite.moran+1
  }
  
  return(H)
}

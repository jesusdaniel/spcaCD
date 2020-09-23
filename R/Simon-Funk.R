
library('Matrix')
library('irlba')




iter.SVD.core <- function(A,K,tol=1e-5,max.iter=100,sparse=TRUE,init=NULL,verbose=FALSE,tau=0,fast=FALSE,p.sample=1,kappa=NULL){
  #browser()
  if(sparse) A <- Matrix(A,sparse=TRUE)
  avg.p <- mean(as.numeric(A),na.rm=TRUE)
  if(is.null(kappa))kappa <- 2/avg.p
  cap <- 1#kappa*avg.p
  if(cap>1-tau) cap <- 1-tau
  if(fast){
    if(verbose) print("Matrix completion with fast approximation!")
    A[which(is.na(A))] <- 0
    A <- A/p.sample
    #svd.new <- svd(A,nu=K,nv=K)
    svd.new <- irlba(A,nu=K,nv=K)
    if(K==1){ A.new <- svd.new$d[1]*matrix(svd.new$u,ncol=1)%*%t(matrix(svd.new$v,ncol=1))}else{
      A.new <- svd.new$u%*%(t(svd.new$v)*svd.new$d[1:K])}
    A.new[A.new < 0+tau] <- 0+tau
    A.new[A.new >cap] <- cap
    return(list(iter=NA,SVD=svd.new,A=A.new,err.seq=NA))
  }
  #if(sparse) A <- Matrix(A,sparse=TRUE)
  Omega <- which(is.na(A))
  A.col.means <- colMeans(A,na.rm=TRUE)
  #avg.p <- mean(as.numeric(A),na.rm=TRUE)
  A.impute <- A
  n <- nrow(A)
  p <- ncol(A)
  if(is.null(init)){
    A.impute[Omega] <- runif(n=length(Omega))
    init.SVD <- irlba(A.impute,nu=K,nv=K)
  }else{
    init.SVD <- init
  }
  #print(init.SVD$u)
  ## init.SVD <- irlba(A,nu=K,nv=K) ## if you are working on large problems
  if(K==1){U.old <- V.old <- matrix(0,n,1)}else{
    if(K==2){
      U.old <- matrix(init.SVD$u[,1:(K-1)],ncol=K-1)*(sqrt(init.SVD$d[1:(K-1)]))
      V.old <- matrix(init.SVD$v[,1:(K-1)],ncol=K-1)*(sqrt(init.SVD$d[1:(K-1)]))
    }else{
      #print(init.SVD$u)
      U.old <- matrix(init.SVD$u[,1:(K-1)],ncol=K-1)%*%diag(sqrt(init.SVD$d[1:(K-1)]))
      V.old <- matrix(init.SVD$v[,1:(K-1)],ncol=K-1)%*%diag(sqrt(init.SVD$d[1:(K-1)]))
    }
  }
  A.old <- U.old %*% t(V.old)
  R <- A - A.old
  
  R[Omega] <- 0
  if(verbose) print(norm(R))
  A.impute <- A
  A.impute[Omega] <- A.old[Omega]
  ### begin iteration
  flag <- 0
  iter <- 0
  err.seq <- norm(R,"F")
  shrink <- 0
  while((iter < max.iter) && (flag != 1)){
    #print(iter)
    iter <- iter + 1
    svd.new <- irlba(A.impute,nu=K,nv=K)
    if(K==1){ A.new <- svd.new$d[1]*matrix(svd.new$u,ncol=1)%*%t(matrix(svd.new$v,ncol=1))}else{
      A.new <- svd.new$u%*%diag(svd.new$d)%*%t(svd.new$v)}
    A.new[A.new < 0+tau] <- 0+tau
    A.new[A.new >cap] <- cap
    A.impute[Omega] <- A.new[Omega]
    A.old <- A.new
    R <- A.impute - A.new
    err <- norm(R,"F")
    if(verbose) print(err)
    err.seq <- c(err.seq,err)
    if(abs(err.seq[iter+1]-err.seq[iter])<tol*err.seq[iter]) flag <- 1
  }
  #print(iter)
  return(list(iter=iter,SVD=svd.new,A=A.new,err.seq=err.seq))
  #return(list(iter=iter,SVD=svd.new,A=A.impute,err.seq=err.seq))
}

trapezoid <- function(x, y) {
  # numerical integral of fun from a to b
  # using the trapezoid rule with n subdivisions
  # assume a < b and n is a positive integer
  h <- diff(x)
  n <- length(y) - 1
  s <- sum(y[1:n] *h + (y[2:(n+1)] - y[1:n]) *h/2)
  return(s)
}



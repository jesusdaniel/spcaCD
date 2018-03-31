# Likelihood

sbm_loglikelihood <- function(A,Z) {
  if(sum(apply(Z,1,sum))!=nrow(Z))
    return(NA)
  B_sum = crossprod(Z,A)%*%Z
  com_sizes = diag(t(crossprod(Z)))
  B = B_sum / (tcrossprod(com_sizes) - diag(com_sizes))
  W = Z%*%tcrossprod(B,Z)
  Loglik = A*log(W)+(1-A)*log(1-W)
  return(sum(Loglik[upper.tri(Loglik)]))
}

W_loglikelihood <- function(A,W) {
  W=ifelse(W<=0,0.00001,W)
  Loglik = A*log(W)+(1-A)*log(1-W)
  #where = (W>0 & W<1)
  #return(sum(Loglik[upper.tri(Loglik)*where])/(sum(upper.tri(Loglik)*where)/2))
  return(sum(Loglik[upper.tri(Loglik)]))
}

# Get the min squares B
B_min_squares <- function(A,Z) {
  K = ncol(Z)
  n=ncol(A)
  B = array(0,dim = c(K,K))
  B_indexes = which(!lower.tri(B),arr.ind = T)
  XX = apply(B_indexes, 1,function(x1) apply(B_indexes,1, 
                                             function(x2) {
                                               b=(Z[,x1[1]]%*%t(Z[,x1[2]]))*(Z[,x2[1]]%*%t(Z[,x2[2]]))
                                               diag(b) = 0
                                               sum(b)/(1+1*(x1[1]==x2[1]&x1[2]==x2[2] & x1[1]==x1[2]))
                                             }))
  
  XY = apply(B_indexes,1,function(x1) {
    b=(Z[,x1[1]]%*%t(Z[,x1[2]]))*A
    diag(b) = 0
    sum(b)/(1+1*(x1[1]==x1[2]))})

  B[!lower.tri(B)] = solve(XX,XY)
  B[lower.tri(B)] = t(B)[lower.tri(B)]
  W = Z%*%B%*%t(Z)
  loglik = W_loglikelihood(A,W)
  return(list(B=B,W =W,loglik = loglik))
}




# Solve the problem
# minimize_{P,B} -logklihelihood(P)
# Subject to P = ZBZ'
B_fit <- function(A, Z, mu = 10) {
  # Initial B
  n = ncol(A)
  Bk = array(sum(A)/n^2,dim = c(ncol(Z),ncol(Z)))
  Uk= array(0, dim = c(nrow(Z),nrow(Z)))
  B_indexes = which(!lower.tri(Bk),arr.ind = T)
  indexes = which(upper.tri(A),arr.ind= T)
  
  Y = A[upper.tri(A)]

  X = t(apply(indexes,1,function(uv)     apply(B_indexes, 1,function(ij) if(ij[1]==ij[2]){
    Z[uv[1],ij[1]]*Z[uv[2],ij[2]] } else{ Z[uv[1],ij[1]]*Z[uv[2],ij[2]] + Z[uv[1],ij[2]]*Z[uv[2],ij[1]]})))
  
  crossprod(X)
  crossprod(X,Y)
  
  B_indexes = which(!lower.tri(Bk),arr.ind = T)
  XX = apply(B_indexes, 1,function(x1) apply(B_indexes,1, 
                                             function(x2) {
                                               b=(Z[,x1[1]]%*%t(Z[,x1[2]]))*(Z[,x2[1]]%*%t(Z[,x2[2]]))
                                               diag(b) = 0
                                               sum(b)/(1+1*(x1[1]==x2[1]&x1[2]==x2[2] & x1[1]==x1[2]))
                                             }))
  
  XY = apply(B_indexes,1,function(x1) {
    b=(Z[,x1[1]]%*%t(Z[,x1[2]]))*A
    diag(b) = 0
    sum(b)/(1+1*(x1[1]==x1[2]))})
  XX
  solve(XX,XY)
  
  t(Z)%*%A%*%Z
  t(Z)%*%Z
  
  V = Z%*%t(Z)
  diag(V)=0
  
  b = ((Z[,1]%*%t(Z[,2]))*(Z[,2]%*%t(Z[,2])))
  
  diag(b) = 0
  sum(b)
  
  sum(Z[,2]*Z[,2]*Z[,1]*Z[,1])
  
  V = Z[,2]%*%t(Z[,1])
  
  diag(V) = 0
  sum(V^2)/2
  
  det(t(Z)%*%Z)
  
  
  for( i in 1:10) {
  # Pk step
  H = Z%*% tcrossprod(Bk,Z)-Uk
  Pk = (H+A-1)/2 + sqrt((H+A-1)^2/4 + 1/mu)
  #levelplot(Pk)
  #Bk step
  M = Pk + Uk
  Y = M[upper.tri(M)]
  

  Bk_lm = lm(Y~X-1)
  Bk[!lower.tri(Bk)] = Bk_lm$coefficients
  Bk[!upper.tri(Bk)] = Bk_lm$coefficients
  Uk = Uk + Pk - Z%*%tcrossprod(Bk,Z)
  }
  Bk
  mean((A[1:250,1:250])[upper.tri((A[1:250,1:250]))])
  + 1/250
  mean((A[1:239]))
}

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
  W= pmin(pmax(W, 1e-3), 1-1e-3)
  Loglik = A*log(W)+(1-A)*log(1-W)
  #where = (W>0 & W<1)
  #return(sum(Loglik[upper.tri(Loglik)*where])/(sum(upper.tri(Loglik)*where)/2))
  return(sum(Loglik[upper.tri(Loglik)]))
}

get_W_estimate <- function(A, Z) {
  W = as.matrix((Z%*%(ginv(as.matrix(Z))%*% (A%*% t(ginv(as.matrix(Z)))))) %*% t(Z))
  W = ifelse(W>1,0.999,W); W = ifelse(W<0,0.001,W)
}

loglikelihood <- function(A, Z) {
  What <- get_W_estimate(A, Z)
  W_loglikelihood(A, What)
}

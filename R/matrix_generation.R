# Matrix generation

generate_overlapping <- function(rho = 0.1, d = 20,   n = 2000) {
  require(Matrix)
  K = 3
  pi1 = 0.3#0.25
  pi2 = 0.03#0.05
  pi3 = 0.01#0.1
  B = array(rho,dim = c(K,K))
  diag(B) = 1
  Z = Matrix(0,nrow = n, ncol=K)
  Z[1:(n*pi1),1] = 1
  Z[(n*pi1+1):(2*n*pi1),2] = 1
  Z[(2*n*pi1+1):(3*n*pi1),3] = 1
  Z[(3*n*pi1 + 1):(3*n*pi1 + n*pi2),c(1,2)] = 1/sqrt(2)
  Z[(3*n*pi1 + n*pi2 + 1):(3*n*pi1 + 2*n*pi2),c(2,3)] = 1/sqrt(2)
  Z[(3*n*pi1 + 2*n*pi2+ 1):(3*n*pi1 + 3*n*pi2),c(1,3)] = 1/sqrt(2)
  Z[(3*n*pi1 + 3*n*pi2+1):n,] = 1/sqrt(3)
  A = Z%*%B%*%t(Z)
  alpha = sum(apply(A,2,sum)-diag(A))/(n*d)
  W = A/alpha
  Aobs = apply(W,MARGIN = c(1,2),function(u) rbinom(1,1,prob = u))
  
  Aobs[lower.tri(Aobs)] = 0
  Aobs =  Aobs + t(Aobs)
  diag(Aobs) = 0
  return(list(Aobs= Aobs, membership = Z,B = B, W = W, alpha = alpha))
}

laplacian <- function(A, tau = 0) {
  Atau = A+tau
  degree = apply(Atau,1,sum)
  diag(1/sqrt(degree))%*%Atau%*%diag(1/sqrt(degree))
}


generate_overlapping_hubs <- function(rho = 0.1, d = 20,   n = 500, p_hub = 0.2, hub_value = sqrt(20)) {
  require(Matrix)
  K = 3
  pi1 = 0.3#0.25
  pi2 = 0.03#0.05
  pi3 = 0.01#0.1
  Theta = rbinom(n,1,p_hub)*(hub_value-1)+1
  Theta = Theta/sum(Theta)*n
  B = array(rho,dim = c(K,K))
  diag(B) = 1
  Z = Matrix(0,nrow = n, ncol=K)
  Z[1:(n*pi1),1] = 1
  Z[(n*pi1+1):(2*n*pi1),2] = 1
  Z[(2*n*pi1+1):(3*n*pi1),3] = 1
  Z[(3*n*pi1 + 1):(3*n*pi1 + n*pi2),c(1,2)] = 1/sqrt(2)
  Z[(3*n*pi1 + n*pi2 + 1):(3*n*pi1 + 2*n*pi2),c(2,3)] = 1/sqrt(2)
  Z[(3*n*pi1 + 2*n*pi2+ 1):(3*n*pi1 + 3*n*pi2),c(1,3)] = 1/sqrt(2)
  Z[(3*n*pi1 + 3*n*pi2+1):n,] = 1/sqrt(3)
  A = diag(Theta)%*%Z%*%B%*%t(Z)%*%diag(Theta)
  diag(A) = 0
  alpha = sum(apply(A,2,sum))/(n*d)
  W = A/alpha
  Aobs = apply(W,MARGIN = c(1,2),function(u) rbinom(1,1,prob = u))
  
  Aobs[lower.tri(Aobs)] = 0
  Aobs =  Aobs + t(Aobs)
  diag(Aobs) = 0
  return(list(Aobs= Aobs, membership = Z,B = B, W = W))
}

generate_occam <- function(Z, B, alpha = 1, theta = NULL) {
  K = ncol(Z)
  N = nrow(Z)
  thetaZ <- Z
  if(!is.null(theta))
    thetaZ = diag(theta) %*% Z
  W <- thetaZ %*% (alpha * B) %*% t(thetaZ)
  W <- pmax(pmin(W, 1), 0)
  Aup <- sapply(W[upper.tri(W)], function(x) rbinom(1, 1, prob = x))
  A <- Matrix(0, nrow = N, ncol = N)
  A[upper.tri(A)] <- Aup
  A =  A + t(A)
  return(list(A = A, W = W))
}


generate_overlapping_k3 <- function(rho = 0.1, d = 20, npure = 450, 
                                    noverlap2 = 45, noverlap3 = 5, phubs = 0, hubsize= 5) {
  require(Matrix)
  K = 3
  n = N = npure + noverlap2 + noverlap3
  pi1 = npure/(N*K)
  pi2 = noverlap2/(N*K)
  pi3 = noverlap3/(N)
  B = array(rho,dim = c(K,K))
  diag(B) = 1
  
  Z = Matrix(0, nrow = N, ncol=K)
  
  Z[1:(n*pi1),1] = 1
  Z[(n*pi1+1):(2*n*pi1),2] = 1
  Z[(2*n*pi1+1):(3*n*pi1),3] = 1
  if(noverlap2 > 0) {
    Z[(3*n*pi1 + 1):(3*n*pi1 + n*pi2),c(1,2)] = 1/(2)
    Z[(3*n*pi1 + n*pi2 + 1):(3*n*pi1 + 2*n*pi2),c(2,3)] = 1/(2)
    Z[(3*n*pi1 + 2*n*pi2+ 1):(3*n*pi1 + 3*n*pi2),c(1,3)] = 1/(2) 
  }
  if(noverlap3 > 0) {
    Z[(3*n*pi1 + 3*n*pi2+1):n,] = 1/(3) 
  }
  
  
  
  #
  # generate data
  if(phubs > 0){
    ishub = 1*(runif(N) <= phubs)
    theta = rep(1, N)  + (hubsize-1)*ishub
    Zsum = colSums(diag(theta)%*%Z)
    alpha = (N*d) / (as.double(crossprod(Zsum, B) %*% Zsum)- N)
    Aobs <- generate_occam(Z, B, alpha = alpha, theta = theta)
  }else{
    Zsum = colSums(Z)
    alpha = (N*d) / as.double(crossprod(Zsum, B) %*% Zsum)
    Aobs <- generate_occam(Z, B, alpha = alpha)
  }
  
  
  
  return(list(Aobs= Aobs$A, membership = Z, W = Aobs$W))
}

sbm <- function(p,q, K, nk) {
  Psi = (p-q)*diag(K)  + q
  A = kronecker(Psi,matrix(rep(1,(nk)^2),nrow = nk))
  Aobs = apply(A,MARGIN = c(1,2),function(u) rbinom(1,1,prob = u))  
  Aobs[lower.tri(Aobs)] = 0
  Aobs =  Aobs + t(Aobs)
  diag(Aobs) = 0
  return(list(Aobs= Aobs, membership = kronecker(1:K,rep(1,nk)), W = A))
}

dc_erdosrenyi <- function(p,theta,n) {
  A = p*matrix(rep(1,(n)^2),nrow = n)
  A = diag(theta)%*%A%*%diag(theta)
  Aobs = apply(A,MARGIN = c(1,2),function(u) rbinom(1,1,prob = min(u,0.99)))  
  Aobs[lower.tri(Aobs)] = 0
  Aobs =  Aobs + t(Aobs)
  diag(Aobs) = 0
  return(Aobs)
}

osbm_2 <- function(p,q, n1, n2, n3) {
  Z = array(0, dim = c(n1+n2+n3, K))
  Z[1:n1,1] = 1; Z[(n1+1):(n1+n2),2] = 1; if(n3>0)Z[(n1+n2+1):(n1+n2+n3),] = 1/2; 
  B = (p-q)*diag(K)  + q
  W = Z%*%B%*%t(Z)
  Aobs = apply(W,MARGIN = c(1,2),function(u) rbinom(1,1,prob = u))  
  Aobs[lower.tri(Aobs)] = 0
  Aobs =  Aobs + t(Aobs)
  diag(Aobs) = 0
  return(list(Aobs= Aobs, membership = Z, W = W))
}

sbm_differentsize <- function(p,q, K, nks) {
  Psi = (p-q)*diag(K)  + q
  Z = array(0,dim = c(sum(nks),K))
  index = 0
  for(i in 1:K){
    Z[which((1:nrow(Z))>index & (1:nrow(Z))<=index+nks[i]),i]=1
    index = index + nks[i]
  }
  A = Z%*%Psi%*%t(Z)
  Aobs = apply(A,MARGIN = c(1,2),function(u) rbinom(1,1,prob = u))  
  Aobs[lower.tri(Aobs)] = 0
  Aobs =  Aobs + t(Aobs)
  diag(Aobs) = 0
  return(list(Aobs= Aobs, membership = Z))
}


overlapping_2k <- function(rho = 0.1, d = 20,   n = 100, overlapping_nodes = 1) {
  require(Matrix)
  K = 2
  B = array(rho,dim = c(K,K))
  diag(B) = 1
  Z = Matrix(0,nrow = n+overlapping_nodes, ncol=K)
  Z[1:n/2,1] = 1
  Z[(n/2+1):n,2] = 1
  Z[(n+1):(n+overlapping_nodes),] = 1/2
  A = Z%*%B%*%t(Z)
  alpha = sum(apply(A,2,sum))/(n*d)
  W = A/alpha
  Aobs = apply(W,MARGIN = c(1,2),function(u) rbinom(1,1,prob = u))
  
  Aobs[lower.tri(Aobs)] = 0
  Aobs =  Aobs + t(Aobs)
  diag(Aobs) = 0
  return(list(Aobs= Aobs, membership = Z,B = B, W = W))
}

sbm_cai  <- function() {
  p = 0.17; q = 0.11; K = 2; nk = 500
  Psi = (p-q)*diag(K)  + q
  A = kronecker(Psi,matrix(rep(1,(nk)^2),nrow = nk))
  W = array(0.7, dim = c(30,30))
  beta = runif(1000)^2
  Z = beta%*%t(rep(1,30))
  P = cbind(A,Z)
  Z = cbind(t(Z),W)
  P = rbind(P,Z)
  Aobs = apply(P,MARGIN = c(1,2),function(u) rbinom(1,1,prob = u))  
  Aobs[lower.tri(Aobs)] = 0
  Aobs =  Aobs + t(Aobs)
  diag(Aobs) = 0
  return(list(Aobs= Aobs))
}

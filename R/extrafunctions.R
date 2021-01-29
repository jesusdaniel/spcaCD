tolerance_function = function(u,v) {
  #nonzeros = sum(u -  v %*% diag(sign(diag(crossprod(u,v))))!=0)
  #if(nonzeros==0)
  #  return(0)
  #norm(u -  v %*% diag(sign(diag(crossprod(u,v)))))/sqrt(nonzeros)
  norm(u-v)/sum(u!=0)
}

leading_eigenspace <- function(A, K) {
  eig <- eigen(A)
  order = order(abs(eig$values), decreasing = T)
  return(list(V = eig$vectors[,order[1:K]], lambdas = eig$values[order[1:K]]))
}

initialization = function(N, K) {
  initial = array(0, dim = c(N,K))
  flag = FALSE
  while(!flag) {
    indexes = floor(1+runif(min = 0,max = K, n = N))
    for(i in 1:K)
      initial[which(indexes==i),i] = 1
    sums_cols = colSums(initial)
    flag = min(sums_cols)/max(sums_cols) >0.3
  }
  return(initial)
}

initialization_kmeans <- function(A, K) {
  #browser()
  eig <- eigen(A)
  #Vs = eig$vectors[, order(abs(eig$values), decreasing = T)[1:K]]
  Vs = eig$vectors[, 1:K]
  if(sum(Vs[,1] ==0)==0) {
    Rs <- Vs[,2:K] / Vs[,1]
    Rsthresh <- pmax(Rs, -log(ncol(A)))
    Rsthresh <- pmin(Rsthresh, log(ncol(A)))
    km <- kmeans(Rsthresh, K, nstart = 10)
  }else {
    km <- kmeans(Vs, K, nstart = 10)
  }
  Zinit <- matrix(0, ncol(A), K)
  for(i in 1:K) {
    Zinit[which(km$cluster==i),i] = 1
  }
  Zinit
}

# reconstruction error of ||A - P_Z(A)||^2
Frobenius_error <- function(Z, A) {
  sum((ginv(crossprod(as.matrix(Z))) %*% (crossprod(Z, A) %*% Z) )^2)
}
PCA_trace = function(Z, A) {
  return(sum(diag(ginv((crossprod(as.matrix(Z))))) %*% ((t(Z)%*%A)%*%Z)))
}

PCA_trace_degree = function(Z, A) {
  V = abs(eigen(A)$vectors[,1])
  V = V/sum(V)
  return(sum(diag(solve(crossprod(diag(V)%*%Z)) %*% ((t(diag(V)%*%Z)%*%A)%*%(diag(V)%*%Z)))))
}
fit_Zest_l1 = function(Z) {
  t(apply(Z,1,function(z)  z/(sum(abs(z))+ 1*(sum(abs(z))==0))))
}

laplacian <- function(A, tau = 0) {
  Atau = A+tau
  degree = apply(Atau,1,sum)
  diag(1/sqrt(degree))%*%Atau%*%diag(1/sqrt(degree))
}


fit_Zest_l2 = function(Z) {
  t(apply(Z,1,function(z)  z/(sqrt(sum(abs(z)^2))+ 1*(sum(abs(z))==0))))
}


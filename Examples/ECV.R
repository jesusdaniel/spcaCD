## Cross validation by matrix completion
sample_nodepairs <- function(N, p) {
  selected = matrix(0, N, N)
  selected[upper.tri(selected)] = rbinom(N*(N-1)/2, size = 1, p)
  selected <- selected + t(selected)
  selected
}

svd_impute <- function(mat, K){
  eigen = eigen(mat)
  ord = order(abs(eigen$values), decreasing = T)
  return(eigen$vectors[,ord[1:K]] %*% 
           diag(eigen$values[ord[1:K]], nrow = K) %*% t(eigen$vectors[,ord[1:K]]))
}

simple_impute <- function(A, Omega){
  diag(Omega) = 0
  phat = sum(A*Omega) / sum(Omega)
  diag(Omega) = 1
  Aimp = A*Omega + phat*(1-Omega)
  Aimp
}

loss_function <- function(A, Omega, Ahat) {
  diag(Omega) = 1
  #Ahat = pmin(pmax(0, Ahat),1)
  sum( abs((A-Ahat)*(1-Omega))/sum(1-Omega))
}

lambda = 0.5
lambdapath = seq(0, 0.7, 0.05)

results <- sapply(lambdapath, function(lambda) {
  M = 10
  cvloss = rep(0, M)
  sparsity = rep(0, M)
  for(i in 1:M) {
    p = 0.9
    N = ncol(A)
    Omega = sample_nodepairs(ncol(A), p)
    Mhat = simple_impute(A, Omega)#svd_impute(Omega*A/p, K)
    #Mhat = svd_impute(Omega*A/p, K)
    Z0 <- Ztrue#initialization(N, K)
    sol_spca <- spca_iterative_DC(Mhat, K, lambda = lambda, lnorm = 2,
                                  Z0 = Z0, thresh = soft_thresh,
                                  num_iter = 50, TOL = 1e-06, save_path = T)
    sparsity[i] = sum(abs(sol_spca$Z)!=0)
    #sol_spca_A <- spca_iterative_DC(A, K, lambda = lambda, lnorm = 2,
    #                              Z0 = Z0, thresh = soft_thresh,
    #                              num_iter = 50, TOL = 1e-06, save_path = T)
    #distance_subs(eigen(Mhat)$vectors[,1:K], sol_spca$Z)
    
    #Zhat = eigen(Mhat)$vectors[,1:K]
    Zhat = sol_spca$Z
    #Phat = Zhat %*% solve(crossprod(Zhat)) %*% t(Mhat %*% (Zhat))
    Phat2 = Zhat%*% solve(crossprod(Zhat)) %*% (crossprod(Zhat, Mhat) %*% Zhat) %*%
      solve(crossprod(Zhat)) %*% t(Zhat)
    
    #cvloss[i] = exNiV(Zhat>0, sol_spca_A$Z>0)#loss_function(A, Omega, Phat2)
    cvloss[i] = loss_function(A, Omega, Phat2)
  } 
  cvloss
  #sparsity
})


rho = 0.4
d = 20
K = 3
A <- generate_overlapping(n = 500, rho = rho, d = d)
Ztrue <- A$membership
A <- A$Aobs
N = nrow(A)
Z0 <- initialization(N, K)
colSums(Ztrue)


M = 10
lambdapath = seq(0, 0.99, 0.03)

results2 <- sapply(1:M, function(i) {
  p = 0.9
  N = ncol(A)
  Omega = sample_nodepairs(ncol(A), p)
  Mhat = simple_impute(A, Omega)
  Mhat = svd_impute(Omega*A/p, K)
  cvloss = rep(0, length(lambdapath))
  sparsity = rep(0, length(lambdapath))
  j=1
  Z0 <- Ztrue
  Z0 <- initialization_kmeans(A, K)
  for(lambda in lambdapath) {
    sol_spca <- spca_iterative_DC(Mhat, K, lambda = lambda, lnorm = 2,
                                  Z0 = Z0, thresh = hard_thresh,
                                  num_iter = 50, TOL = 1e-06, save_path = T)
    sparsity[j] = sum(abs(sol_spca$Z)!=0)
    #sol_spca_A <- spca_iterative_DC(A, K, lambda = lambda, lnorm = 2,
    #                              Z0 = Z0, thresh = soft_thresh,
    #                              num_iter = 50, TOL = 1e-06, save_path = T)
    #distance_subs(eigen(Mhat)$vectors[,1:K], sol_spca$Z)
    
    #Zhat = eigen(Mhat)$vectors[,1:K]
    Zhat = sol_spca$Z
    #Phat = Zhat %*% solve(crossprod(Zhat)) %*% t(Mhat %*% (Zhat))
    Phat2 = Zhat%*% ginv(crossprod(Zhat)) %*% (crossprod(Zhat, Mhat) %*% Zhat) %*%
      ginv(crossprod(Zhat)) %*% t(Zhat)
    
    #cvloss[i] = exNiV(Zhat>0, sol_spca_A$Z>0)#loss_function(A, Omega, Phat2)
    cvloss[j] = loss_function(A, Omega, Phat2)
    j=j+1
  } 
  cvloss
  #plot(cvloss)
  #sparsity
})

plot(apply(results2,1,mean))
points(apply(results2,1,mean) + apply(results2,1,sd)/sqrt(M), col = "red")

lambda = lambdapath[which.min(rowMeans(results2))]
lambda = 0.78
Z0 <- initialization_kmeans(A, K)
Z0 = sol_occam
sol_spca <- spca_iterative_DC(A, K, lambda = lambda, lnorm = 2,
                              Z0 = Z0, thresh = hard_thresh,
                              num_iter = 100, TOL = 1e-06, save_path = T)
distance_subs(eigen(A)$vectors[,1:K], sol_spca$Z)
plot(eigen(A)$vectors[,1])
plot(eigen(A)$vectors[,2])
plot(eigen(A)$values)
exNiV(sol_spca$Z>0, Ztrue>0)
crossprod(binarize(sol_spca$Z), Z)
crossprod((sol_spca$Z)!=0, Z)
crossprod(Z0, Z)
apply(sol_)

plot(sol_spca$Z[,1])
plot(sol_spca$Z[,2])
plot(sol_spca$Z[,3])

plot(Z0[,1])
plot(Z0[,2])
plot(Z0[,3])
which(apply(sol_spca$Z!=0,1, sum)>=4)
sol_spca$Z[73, ]
Z[73,]
Z[,]
names[which(apply(sol_spca$Z!=0,1, sum)>=3)]
partidos_sub[which(apply(sol_spca$Z!=0,1, sum)>=3)]

sum(sol_spca$Z==0)
sum(apply(sol_spca$Z,1,sum)==0)
min(sol_spca$Z)
sum(sol_spca$Z<0)
sol_spca$iterations
plot(sol_spca$tol_path)

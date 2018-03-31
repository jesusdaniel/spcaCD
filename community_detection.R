# Fits overlapping community detection.
# Arguments:
# A = Adjacency matrix
# K = number of communities
# lambda = penalty parameter
# epsilon = everything below this value (in abs. value) will be considered as zero
# iterations = number of iteration of SPCA
# Returns:
# Z = assignments matrix
overlapping_community_detection <- function(A, K, lambda, DC=FALSE, iterations = 100, TOL= 1e-06,
                                            initial_val = NULL, save_path = F,lnorm_sol=1,
                                            epsilon = 1e-10, num_trials = 1) {
 # browser()
  if(is.null(initial_val)) {
    initial_val = initialization(ncol(A), K)
  }
  trials = 1
  if(!DC) {
    sparsePCA = spca_iterative(A, K, lambda, lnorm_sol,
                               initial_val, iterations,TOL = TOL, save_path = save_path)   
  }else{
    sparsePCA = spca_iterative_DC(A, K, lambda, lnorm_sol,
                                  initial_val, iterations,TOL = TOL^2, save_path = save_path)  
  }
  trials = trials + 1
  best_spca  = sparsePCA
  if(num_trials > 1) {
    best_F = PCA_trace(sparsePCA$Z, A)
    while(trials <= num_trials) {
      initial_val = initialization(ncol(A), K)
      if(!DC) {
        sparsePCA = spca_iterative(A, K, lambda, lnorm_sol,
                                   initial_val, iterations,TOL = TOL, save_path = save_path)   
      }else{
        sparsePCA = spca_iterative_DC(A, K, lambda, lnorm_sol,
                                      initial_val, iterations,TOL = TOL^2, save_path = save_path)  
      }
      if(PCA_trace(sparsePCA$Z, A)> best_F) {
        best_F = PCA_trace(sparsePCA$Z, A)
        best_spca = sparsePCA
      }
      trials = trials +1
    } 
  }
  Zraw = best_spca$Z
  Z = abs(best_spca$Z)*(abs(best_spca$Z)> epsilon)
  theta = apply(Z,1,function(x) (sum(abs(x)^lnorm_sol))^(1/lnorm_sol))
  theta = ifelse(theta==0,1,theta)
  Z = Z/theta
  return(list(Zest = Z, Zraw = Zraw, theta = theta,
              membership = 1*(Z!=0), iterations = best_spca$iterations, 
              tol_path = best_spca$tol_path,  Z_path = best_spca$Z_path))
}

PCA_trace = function(Z, A) {
  return(sum(diag(ginv(crossprod(Z)) %*% ((t(Z)%*%A)%*%Z))))
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
  eig <- eigen(A)
  Vs = eig$vectors[, order(abs(eig$values), decreasing = T)[1:K]]
  Rs <- Vs[,2:K] / Vs[,1]
  Zinit <- matrix(0, ncol(A), K)
  km <- kmeans(Rs, K)
  for(i in 1:K) {
    Zinit[which(km$cluster==i),i] = 1
  }
  Zinit
}

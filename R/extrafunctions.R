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



# 
# #### BIC
# get_W_estimate <- function(A, Z) {
#   W = as.matrix((Zest%*%(ginv(Zest)%*% (A%*% t(ginv(Zest))))) %*% t(Zest))
#   W = ifelse(W>1,0.999,W); W = ifelse(W<0,0.001,W)
# }
#       # BIC
#     W = as.matrix((Zest%*%(ginv(Zest)%*% (A%*% t(ginv(Zest))))) %*% t(Zest))
#     W = ifelse(W>1,0.999,W); W = ifelse(W<0,0.001,W)
#     logliks_Zest = c(logliks_Zest, W_loglikelihood(A, W))
#     # Degree corrected BIC
#     #degrees = apply(A,1,sum)
#     #thetas = sqrt(degrees/apply(Zest,1,function(x) sum(x)+1*sum(x==0)))
#     #thetas = thetas/sum(thetas)
#     #A2 = t((1/thetas)*A)*(1/thetas)
#     if(DC){
#       Zest = u$Zraw
#       W = as.matrix(Zest%*%(ginv(Zest)%*% (A%*% t(ginv(Zest))))) %*% t(Zest)
#       W = ifelse(W>1,0.999,W); W = ifelse(W<0,0.001,W)
#       logliks_Zest_DC = c(logliks_Zest_DC, W_loglikelihood(A, W)) 
#     }
#     #logliks_Zest_DC = c(logliks_Zest_DC, W_loglikelihood(A, W))
#     i = i+1
#   }
#   
#   #BIC uncorrected
#   logliks_Zest = ifelse(is.na(logliks_Zest),-Inf,logliks_Zest)
#   logliks_Zest = ifelse(sparsity==0,-Inf,logliks_Zest)
#   BIC_vals = sapply(1:length(logliks_Zest), 
#                     function(i) -2*logliks_Zest[i] + (N*(K-1)-sparsity[i])*log(N))
#   Z_hat_SPCA = path[[which.min(BIC_vals)]]$Zest
#   Z_hat_SPCA_DC = NULL
#   best_accuracy_DC = NA
#   best_exNVI_DC = NA
#   #BIC DC
#   if(DC) {
#     logliks_Zest_DC = ifelse(is.na(logliks_Zest_DC),-Inf,logliks_Zest_DC)
#     logliks_Zest_DC = ifelse(sparsity==0,-Inf,logliks_Zest_DC)
#     #BIC
#     BIC_vals_DC = sapply(1:length(logliks_Zest_DC), 
#                          function(i) -2*logliks_Zest_DC[i] + (N*(K-1)-sparsity[i])*log(N))
#     Z_hat_SPCA_DC = path[[which.min(BIC_vals_DC)]]$Zest
#     if(!is.null(Ztrue)) best_accuracy_DC = accuracy_estimated[which.min(BIC_vals_DC)]
#     if(!is.null(Ztrue))  best_exNVI_DC =  exNVIs_estimated[which.min(BIC_vals_DC)]
#   }
#   return(list( logliks_Zest= logliks_Zest, BIC_vals = BIC_vals,
#                path = path,
#                accuracy_estimated = accuracy_estimated,
#                exNVIs_estimated = exNVIs_estimated,
#                sparsity = sparsity,
#                total_iterations = total_iterations,
#                Z_hat_SPCA = Z_hat_SPCA,
#                Z_hat_SPCA_DC= Z_hat_SPCA_DC,
#                best_accuracy = ifelse(!is.null(Ztrue), accuracy_estimated[which.min(BIC_vals)], Inf),
#                best_accuracy_DC = ifelse(!is.null(Ztrue), best_accuracy_DC, Inf),
#                best_exNVI =  ifelse(!is.null(Ztrue), exNVIs_estimated[which.min(BIC_vals)], Inf),
#                best_exNVI_DC = ifelse(!is.null(Ztrue), best_exNVI_DC, Inf)))
# }


sparsePCA_comdet_path <- function(A, K, lambdapath, DC = FALSE,
                                  Ztrue = NULL, l_norm = 1, num_trials = 1, maxiter = 300,TOL=1e-04) {
  #browser()
  N = ncol(A)
  accuracy_estimated = c()
  exNVIs_estimated = c()
  sparsity = c()
  total_iterations = c()
  logliks_Zest = c()
  logliks_Zest_DC = c()
  path = list()
  i = 1
  accuracy_estimated  = NULL
  exNVIs_estimated = NULL
  for(l in lambdapath) {
    u = overlapping_community_detection(A = A, K = K, lambda = l,DC=DC, iterations = maxiter, TOL = TOL,
                                        initial_val = NULL, save_path = F,lnorm_sol=l_norm, num_trials = num_trials)
    path[[i]] = u
    if(!is.null(Ztrue)) accuracy_estimated = c(accuracy_estimated, accuracy(binarize(u$Zest),binarize(Ztrue)))
    if(!is.null(Ztrue))  exNVIs_estimated = c(exNVIs_estimated,exNiV(u$membership,binarize(Ztrue)) )
    sparsity = c(sparsity, sum(u$membership==0))
    total_iterations = c(total_iterations, u$iterations)
    Zest =u$Zest
    # BIC
    W = as.matrix((Zest%*%(ginv(Zest)%*% (A%*% t(ginv(Zest))))) %*% t(Zest))
    W = ifelse(W>1,0.999,W); W = ifelse(W<0,0.001,W)
    logliks_Zest = c(logliks_Zest, W_loglikelihood(A, W))
    # Degree corrected BIC
    #degrees = apply(A,1,sum)
    #thetas = sqrt(degrees/apply(Zest,1,function(x) sum(x)+1*sum(x==0)))
    #thetas = thetas/sum(thetas)
    #A2 = t((1/thetas)*A)*(1/thetas)
    if(DC){
      Zest = u$Zraw
      W = as.matrix(Zest%*%(ginv(Zest)%*% (A%*% t(ginv(Zest))))) %*% t(Zest)
      W = ifelse(W>1,0.999,W); W = ifelse(W<0,0.001,W)
      logliks_Zest_DC = c(logliks_Zest_DC, W_loglikelihood(A, W)) 
    }
    #logliks_Zest_DC = c(logliks_Zest_DC, W_loglikelihood(A, W))
    i = i+1
  }
  
  #BIC uncorrected
  logliks_Zest = ifelse(is.na(logliks_Zest),-Inf,logliks_Zest)
  logliks_Zest = ifelse(sparsity==0,-Inf,logliks_Zest)
  BIC_vals = sapply(1:length(logliks_Zest), 
                    function(i) -2*logliks_Zest[i] + (N*(K-1)-sparsity[i])*log(N))
  Z_hat_SPCA = path[[which.min(BIC_vals)]]$Zest
  Z_hat_SPCA_DC = NULL
  best_accuracy_DC = NA
  best_exNVI_DC = NA
  #BIC DC
  if(DC) {
    logliks_Zest_DC = ifelse(is.na(logliks_Zest_DC),-Inf,logliks_Zest_DC)
    logliks_Zest_DC = ifelse(sparsity==0,-Inf,logliks_Zest_DC)
    #BIC
    BIC_vals_DC = sapply(1:length(logliks_Zest_DC), 
                         function(i) -2*logliks_Zest_DC[i] + (N*(K-1)-sparsity[i])*log(N))
    Z_hat_SPCA_DC = path[[which.min(BIC_vals_DC)]]$Zest
    if(!is.null(Ztrue)) best_accuracy_DC = accuracy_estimated[which.min(BIC_vals_DC)]
    if(!is.null(Ztrue))  best_exNVI_DC =  exNVIs_estimated[which.min(BIC_vals_DC)]
  }
  return(list( logliks_Zest= logliks_Zest, BIC_vals = BIC_vals,
               path = path,
               accuracy_estimated = accuracy_estimated,
               exNVIs_estimated = exNVIs_estimated,
               sparsity = sparsity,
               total_iterations = total_iterations,
               Z_hat_SPCA = Z_hat_SPCA,
               Z_hat_SPCA_DC= Z_hat_SPCA_DC,
               best_accuracy = ifelse(!is.null(Ztrue), accuracy_estimated[which.min(BIC_vals)], Inf),
               best_accuracy_DC = ifelse(!is.null(Ztrue), best_accuracy_DC, Inf),
               best_exNVI =  ifelse(!is.null(Ztrue), exNVIs_estimated[which.min(BIC_vals)], Inf),
               best_exNVI_DC = ifelse(!is.null(Ztrue), best_exNVI_DC, Inf)))
}
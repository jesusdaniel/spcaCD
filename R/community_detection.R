# Fits SPCA to perform community detection.
# Arguments:
# A = Adjacency matrix
# K = number of communities
# lambda = penalty parameter
# method = Type of SPCA method (SPCAeig, SPCA_Z, SPCA_Z_inv, SPCA_DCZ, SPCA_DCZ_inv)
# initial = initial solution. Options: 1) Random, 2) Given, 3) SCORE
# Z0 = if initial = 2, then Z0 is taken as initial solution
# num_trials = If initial is random, number of random starts
# num_iter, TOL, thresh, lnorm = SPCA parameters
# epsilon = Tolerance for thresholding to zero
# Returns:
# Z = assignments matrix
overlapping_community_detection <- function(A, K, lambda = 0.9, method  = "SPCA_DCZ", 
                                            initial = 1, Z0 = NULL, num_trials = 1,
                                            num_iter = 30, TOL= 1e-04, thresh = hard_thresh, lnorm=1,
                                            epsilon = 1e-10) {
  N = ncol(A)
  # Method
  if(method == "SPCAeig") {
    spca <- spca_eigenbasis
  }else{ if(method== "SPCA_Z") {
    spca <- spca_Z
  }else{ if(method == "SPCA_Z_inv") {
    spca <- spca_Z_inv
  }else{ if(method == "SPCA_DCZ") {
    spca <- spca_DCZ
  }else{ if(method == "SPCA_DCZ_inv") {
    spca <- spca_DCZ_inv
  }}}}}
  
  # Initial solution --------------------------------------
  if(initial == 1) {
    Z0 = initialization(N, K)
  }else{ if(initial ==2) {
    if(ncol(Z0) != K | nrow(Z0) != N) {
      stop("The initial value has wrong dimensions.")
    }
  }else{if(initial == 3) {
    Z0 <- initialization_kmeans(A, K)
  } else{ stop("Wrong initial value.")}}}
  
  # SPCA -------------------------------------------------
  trials = 1
  sparsePCA <- spca(A = A, K = K, lambda = lambda, lnorm = lnorm, 
                    Z0 = Z0, num_iter = num_iter, TOL = TOL, save_path = F,
                    thresh = thresh)
  best_spca  = sparsePCA
  
  # More random initial values? -------------------------
  if(num_trials > 1 & initial == 1) {
    # Loss function
    #best_F = Frobenius_error(sparsePCA$Z, A)
    best_F = distance_L2(A, sparsePCA$Z)
    while(trials <= num_trials) {
      Z0 = initialization(ncol(A), K)
      sparsePCA <- spca(A = A, K = K, lambda = lambda, lnorm = lnorm, 
                        Z0 = Z0, num_iter = num_iter, TOL = TOL, save_path = F,
                        thresh = thresh)
      # Compare loss functions
      newF = distance_L2(A, sparsePCA$Z)
      if(newF < best_F) {
        best_F <- newF
        best_spca = sparsePCA
      }
      trials = trials +1
    } 
  }
  Zraw = best_spca$Z
  tot_iterations <- best_spca$iterations
  Z = (best_spca$Z)*(abs(best_spca$Z)> epsilon)
  theta = apply(Z,1,function(x) (sum(abs(x)^lnorm))^(1/lnorm))
  theta = ifelse(theta==0,1,theta)
  Z = abs(Z)/theta
  return(list(Zest = Z, Zraw = Zraw, theta = theta,
              membership = 1*(Z!=0), Z0 = best_spca$Z_path[[1]], method = method,
              tot_iterations = tot_iterations))
}

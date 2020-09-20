############
# Function for calculating a path of solutions for 
# different threshold values in CD
# A = Adjacency matrix
# K = number of kommunities
# method = Type of SPCA method (SPCAeig (Algorithm 1) or SPCA_Z (Algorithm 2))
# lambdapath = set of threshold parameters
# initial = initial solution. Options: 1) Random, 2) User-specified, 3) SCORE
# Z0 = if initial = 2, then Z0 is taken as initial solution
# num_trials = If initial is random, number of random starts
# Z0_previous = If true, start each run with the previous solution on the path of parameters. only works with initial 2 and 3
# p 2
# num_iter, TOL, thresh, lnorm = SPCA parameters
# epsilon = Tolerance for thresholding to zero in community detection
# Returns:
# 
SPCA_CD_path <- function(A, K, method = "SPCA_Z", 
                         lambdapath = c(0.9),
                         initial = 1, Z0 = NULL, num_trials = 1, Z0_previous = FALSE,
                         num_iter = 30, TOL= 1e-04, thresh = hard_thresh, lnorm=1, epsilon = 1e-10) {
  N = ncol(A)

  # Initial solution --------------------------------------
  if(initial == 1) {
    Z0 = NULL
  }else{ if(initial ==2) {
    if(ncol(Z0) != K | nrow(Z0) != N) {
      stop("The initial value has wrong dimensions.")
    }
  }else{if(initial == 3) {
    Z0 <- initialization_kmeans(A, K)
    initial = 2
  } else{ stop("Wrong initial value.")}}}
  
  Z_prev <- Z0
  
  # Iterations -------------------------------------------
  Z_path = list()
  for(i in 1:length(lambdapath)) {
    lambda = lambdapath[i]
    
    if(Z0_previous) {
      Z_start = Z_prev
    }else {
      Z_start = Z0 
    }
    spcaCD <- overlapping_community_detection(A = A, K = K, lambda = lambda, method  = method, 
                                                initial = initial, Z0 = Z_start, num_trials = num_trials,
                                                num_iter = num_iter, TOL= TOL, thresh = thresh, lnorm = lnorm,
                                                epsilon = epsilon) 

    Z_path[[i]] = spcaCD
    Z_prev <- spcaCD$Zraw
  }
  valid_memberships <- sapply(Z_path, function(x) 
    sum(abs(colSums(x$Zraw)) == colSums(abs(x$Zraw)))==K & # all positives or negatives
      sum(colSums(abs(x$Zraw))>0)==K & # nonzero colums
      sum(colSums(x$Zraw==0)>=1) ==K)  # existence of pure nodes
  if(sum(valid_memberships)==0) {
    valid_memberships <- valid_memberships | !valid_memberships
  }
  return(list(Z_path = Z_path, lambdapath = lambdapath, valid_memberships = valid_memberships))
}

SPCA_CD_path.BICselection <- function(A, SPCApath_obj) {
  K <- ncol(SPCApath_obj$Z_path[[1]]$Zraw)
  valid_memberships <- sapply(SPCApath_obj$Z_path, function(x) 
    sum(abs(colSums(x$Zraw)) == colSums(abs(x$Zraw)))==K &
      sum(colSums(abs(x$Zraw))>0)==K)
  if(sum(valid_memberships) == 0)
    valid_memberships <- valid_memberships | !valid_memberships
  lambdas <- SPCApath_obj$lambdapath[valid_memberships]
  BIC <- spca_path_BIC(A, SPCApath_obj$Z_path[valid_memberships]) 
  
  return(list(Z_minBIC= (SPCApath_obj$Z_path[valid_memberships])[[which.min(BIC)]], 
              lambda_minBIC = lambdas[[which.min(BIC)]]))
}
## Different measures of error and fit
spca_path_sparsity <- function(Z_path) {
  sapply(Z_path, function(x) sum(x$membership == 0))
}
spca_path_accuracy <- function(Z_path, Ztrue) {
  sapply(Z_path, function(x) accuracy(x$Zest, Ztrue))
}
spca_path_overlapping_accuracy <- function(Z_path, Ztrue) {
  sapply(Z_path, function(x) overlapping_accuracy(x$membership, Ztrue))
}
spca_path_NVI <- function(Z_path, Ztrue) {
  sapply(Z_path, function(x) exNVI(x$membership, Ztrue))
}
spca_path_subspace_distance_error <- function(Z_path, Ztrue) {
  sapply(Z_path, function(x) distance_subs(x$Zraw, Ztrue))
}
spca_path_subspace_eigenspace_dist <- function(A, Z_path) {
  V = leading_eigenspace(A, ncol(Z_path[[1]]$Zraw))$V
  sapply(Z_path, function(x) distance_subs(x$Zraw, V))
}
spca_path_L2_dist <- function(A, Z_path) {
  sapply(Z_path, function(x) distance_L2(A, x$Zraw))
}
spca_path_loglikehood <- function(A, Z_path) {
  sapply(Z_path, function(x) loglikelihood(A, x$Zraw))
}
spca_path_BIC <- function(A, Z_path) {
  N <- ncol(A)
  K = ncol(Z_path[[1]]$Zraw)
  logliks <- spca_path_loglikehood(A, Z_path)
  sparsity <- spca_path_sparsity(Z_path)
    -2*logliks + 2*(N*(K)-sparsity)*log(N)
}

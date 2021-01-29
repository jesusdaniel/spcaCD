############
# Function for calculating a path of solutions for
# different threshold values in CD
# A = Adjacency matrix
# K = number of kommunities
# method = Type of SPCA method (SPCAeig, SPCAeig2 SPCA_Z or SPCA_DCZ, )
# loss_function for cross-validation: 1 - L1, 2 - L2, 3 - Binomial deviance, 4 (not implemented yet) - AUC
# impute_function for matrix completion. 1-simple, 2-svd, 3 - svd iterative
# lambdapath = set of threshold parameters
# folds, num_folds = parameters for cross-validation
# initial = initial solution. Options: 1) Random, 2) Given, 3) SCORE
# Z0 = if initial = 2, then Z0 is taken as initial solution
# num_trials = If initial is random, number of random starts
# Z0_previous = If true, start each run with the previous solution on the path of parameters. only works with initial 2 and 3
# num_iter, TOL, thresh, lnorm = SPCA parameters
# epsilon = Tolerance for thresholding to zero in community detection
# Returns:
#
cv.SPCA_CD_path <- function(A, K, method = "SPCA_Z",
                            loss_function = 1, impute_function = 2,
                            lambdapath = c(0.9), folds = NULL, num_folds = 10,
                            initial = 1, Z0 = NULL, num_trials = 5, Z0_previous = FALSE,
                            num_iter = 30, TOL= 1e-04, thresh = hard_thresh, lnorm=1, epsilon = 1e-10) {
  N = ncol(A)
  if(loss_function ==1) {
    loss_function = loss_L1
  }else{if(loss_function == 2) {
    loss_function <- loss_L2
  }else{if(loss_function == 3) {
    loss_function <- loss_binomialdeviance
  } else{ stop("Incorrect loss fucntion.")}}}

  if(is.null(folds)) {
    require(cvTools)
    ones <- which(as.matrix(A)==1 & upper.tri(A), arr.ind = T)
    zeros <- which(as.matrix(A)==0 & upper.tri(A), arr.ind = T)
    fold_partition1 = cvFolds(n = nrow(ones), K = num_folds)
    fold_list1 <- lapply(1:fold_partition1$K, function(j)
      ones[fold_partition1$subsets[which(fold_partition1$which==j)],])
    fold_partition0 = cvFolds(n = nrow(zeros), K = num_folds)
    fold_list0 <- lapply(1:fold_partition0$K, function(j)
      zeros[fold_partition0$subsets[which(fold_partition0$which==j)],])
    fold_list <- mapply(rbind, fold_list1, fold_list0, SIMPLIFY = F)
  }
  Omegas <- lapply(fold_list, function(fold) {
    Omega = matrix(0, N, N)
    Omega[fold] = 1
    Omega <- Omega + t(Omega)
    Omega <- 1 - Omega
    Omega
  })
  # Matrix completion
  Mhatlist <- lapply(Omegas, function(Omega) {
    if(impute_function==1)
    {
      simple_impute(A, Omega)
    }else{if(impute_function==2) {
      p = 1-1/num_folds
      Mhat = svd_impute(Omega*A/p, K)
    }else{
      A.new = A
      A.new[Omega==0] = NA
      Mhat <- iter.SVD.core(A.new,K=K)$A
      Mhat <- (Mhat + t(Mhat))/2
      #p = 1-1/num_folds
      #Mhat = svd_iterative_impute(A, Omega, K, p)
    }}
  })
  spca_path <- SPCA_CD_path(A= A, K = K, method = method,
                            lambda = lambdapath, thresh = thresh,
                            initial = initial, Z0 = Z0,  num_trials = num_trials,
                            Z0_previous = Z0_previous, num_iter = num_iter, TOL= TOL,
                            lnorm=lnorm, epsilon = epsilon)
  valid_memberships <- sapply(spca_path$Z_path, function(x)
    sum(abs(colSums(x$Zraw)) == colSums(abs(x$Zraw)))==K & # all positives or negatives
      sum(colSums(abs(x$Zraw))>0)==K & # nonzero colums
      sum(colSums(x$Zraw==0)>=1) ==K)  # existence of pure nodes
  if(sum(valid_memberships)==0) {
    valid_memberships <- valid_memberships | !valid_memberships
  }
  cv.spca <- lapply(Mhatlist, function(Mhat)
    SPCA_CD_path(A = Mhat, K = K, method = method, lambdapath = lambdapath[valid_memberships],
                 initial = initial, Z0 = Z0, num_trials = num_trials, Z0_previous = Z0_previous,
                 num_iter = num_iter, TOL= TOL, thresh = thresh, lnorm=lnorm, epsilon = epsilon))


  #spca_path_NVI(cv.spca[[1]]$Z_path, Ztrue)
  # Evaluate performance
  cv_loss <- sapply(1:length(cv.spca), function(fold) {
    sapply(cv.spca[[fold]]$Z_path, function(sol) {
      Zhat <- as.matrix(sol$Zraw)
      Phat = Zhat%*% ginv(crossprod(Zhat)) %*% (Matrix::crossprod(Zhat, Mhatlist[[fold]]) %*% Zhat) %*%
        ginv(crossprod(Zhat)) %*% t(Zhat)
      loss_function(A, Omegas[[fold]], Phat)
    }) })
  cv_sparsity <- sapply(cv.spca,  function(x) spca_path_sparsity(x$Z_path))
  if(length(cv_loss) > num_folds) {
    minCV = which.min(rowMeans(cv_loss))
  }else{
    minCV=1
  }


  return(list(Z_minCV = (spca_path$Z_path[valid_memberships])[[minCV]],
              lambda_minCV = (lambdapath[valid_memberships])[minCV],
              cv_loss = cv_loss, cv_sparsity = cv_sparsity,
              spca_path = spca_path,
              valid_lambdas = valid_memberships,
              fold_list = fold_list))
}




cv.KSPCA_CD_path <- function(A, Kvals,
                            method = "SPCA_Z",
                            loss_function = 1, impute_function = 2,
                            lambdapath = c(0.9), folds = NULL, num_folds = 10,
                            initial = 1, Z0 = NULL, num_trials = 5, Z0_previous = FALSE,
                            num_iter = 30, TOL= 1e-04, thresh = hard_thresh, lnorm=1,
                            epsilon = 1e-10) {
  N = ncol(A)
  if(loss_function ==1) {
    loss_function = loss_L1
  }else{if(loss_function == 2) {
    loss_function <- loss_L2
  }else{if(loss_function == 3) {
    loss_function <- loss_binomialdeviance
  } else{ stop("Incorrect loss fucntion.")}}}

  if(is.null(folds)) {
    require(cvTools)
    ones <- which(as.matrix(A)==1 & upper.tri(A), arr.ind = T)
    zeros <- which(as.matrix(A)==0 & upper.tri(A), arr.ind = T)
    fold_partition1 = cvFolds(n = nrow(ones), K = num_folds)
    fold_list1 <- lapply(1:fold_partition1$K, function(j)
      ones[fold_partition1$subsets[which(fold_partition1$which==j)],])
    fold_partition0 = cvFolds(n = nrow(zeros), K = num_folds)
    fold_list0 <- lapply(1:fold_partition0$K, function(j)
      zeros[fold_partition0$subsets[which(fold_partition0$which==j)],])
    fold_list <- mapply(rbind, fold_list1, fold_list0, SIMPLIFY = F)
  }
  Omegas <- lapply(fold_list, function(fold) {
    Omega = matrix(0, N, N)
    Omega[fold] = 1
    Omega <- Omega + t(Omega)
    Omega <- 1 - Omega
    Omega
  })

  spca_path_l <- list()
  cv_loss_l <- list()
  cv_sparsity_l <- list()
  res_l <- list()
  minCV_c <- rep(NA, length(Kvals))
  i <- 1
  for(K in Kvals) {
    # Matrix completion
    Mhatlist <- lapply(Omegas, function(Omega) {
      if(impute_function==1)
      {
        #simple_impute(A, Omega)
        A.new = A
        A.new[Omega==0] = NA
        Mhat <- iter.SVD.core(A.new,K=K)$A
        Mhat <- (Mhat + t(Mhat))/2
        Mhat <- A*Omega + Mhat*(1-Omega)
        diag(Mhat) = 0
        Mhat
      }else{if(impute_function==2) {
        p = 1-1/num_folds
        Mhat = svd_impute(Omega*A/p, K)
      }else{
        A.new = A
        A.new[Omega==0] = NA
        Mhat <- iter.SVD.core(A.new,K=K)$A
        Mhat <- (Mhat + t(Mhat))/2
        #p = 1-1/num_folds
        #Mhat = svd_iterative_impute(A, Omega, K, p)
      }}
    })
    spca_path <- SPCA_CD_path(A= A, K = K, method = method,
                              lambda = lambdapath, thresh = thresh,
                              initial = initial, Z0 = Z0,  num_trials = num_trials,
                              Z0_previous = Z0_previous, num_iter = num_iter, TOL= TOL,
                              lnorm=lnorm, epsilon = epsilon)
    spca_path_l[[i]] <- spca_path
    valid_memberships <- sapply(spca_path$Z_path, function(x)
      sum(abs(colSums(x$Zraw)) == colSums(abs(x$Zraw)))==K & # all positives or negatives
        sum(colSums(abs(x$Zraw))>0)==K & # nonzero colums
        #sum(colSums(x$Zraw[which(apply(x$Zraw!=0, 1, function(x)
        #  sum(x)==1)),])==0)==0 & # existence of pure nodes
        sum(2^(1:K) %in% c((x$Zraw != 0) %*% 2^(1:K)))==K &
        !any(duplicated(x$Zraw != 0, MARGIN = 2))) # different columns

    if(sum(valid_memberships)==0) {
      valid_memberships <- 1
    }
    cv.spca <- lapply(Mhatlist, function(Mhat)
      SPCA_CD_path(A = Mhat, K = K, method = method, lambdapath = lambdapath[valid_memberships],
                   initial = initial, Z0 = Z0, num_trials = num_trials, Z0_previous = Z0_previous,
                   num_iter = num_iter, TOL= TOL, thresh = thresh, lnorm=lnorm, epsilon = epsilon))
    if(sum(valid_memberships) == 1) {
        cv_loss <- matrix(sapply(1:length(cv.spca), function(fold) {
          sapply(cv.spca[[fold]]$Z_path, function(sol) {
            Zhat <- as.matrix(sol$Zraw)
            Phat = Zhat%*% ginv(crossprod(Zhat)) %*% (Matrix::crossprod(Zhat, Mhatlist[[fold]]) %*% Zhat) %*%
              ginv(crossprod(Zhat)) %*% t(Zhat)
            loss_function(A, Omegas[[fold]], Phat)
          }) }), nrow = 1)
        cv_sparsity <- matrix(sapply(cv.spca,  function(x) spca_path_sparsity(x$Z_path)), nrow = 1)
      } else {
        cv_loss <- sapply(1:length(cv.spca), function(fold) {
          sapply(cv.spca[[fold]]$Z_path, function(sol) {
            Zhat <- as.matrix(sol$Zraw)
            Phat = Zhat%*% ginv(crossprod(Zhat)) %*% (Matrix::crossprod(Zhat, Mhatlist[[fold]]) %*% Zhat) %*%
              ginv(crossprod(Zhat)) %*% t(Zhat)
            loss_function(A, Omegas[[fold]], Phat)
          }) })
        cv_sparsity <- sapply(cv.spca,  function(x) spca_path_sparsity(x$Z_path))
      }




    if(length(cv_loss) > num_folds) {
      minCV = which.min(rowMeans(cv_loss))
    }else{
      minCV= 1
    }
    cv_loss_l[[i]] <- cv_loss
    spca_path_l[[i]] <- spca_path

    minCV_c[i] <- rowMeans(cv_loss)[minCV]

    res_l[[i]] <- list(Z_minCV = (spca_path$Z_path[valid_memberships])[[minCV]],
         lambda_minCV = (lambdapath[valid_memberships])[minCV],
         cv_loss = cv_loss, cv_sparsity = cv_sparsity,
         spca_path = spca_path,
         valid_lambdas = valid_memberships,
         fold_list = fold_list)

    i <- i + 1
  }
  # Choose K
  best_K <- Kvals[which.min(minCV_c)]


  return(list(best_K = best_K, res_list = res_l))
}

svd_impute <- function(mat, K){
  eigen = eigen(mat)
  ord = order(abs(eigen$values), decreasing = T)
  return(eigen$vectors[,ord[1:K]] %*%
           diag(eigen$values[ord[1:K]], nrow = K) %*% t(eigen$vectors[,ord[1:K]]))
}

svd_iterative_impute <- function(A, Omega, K, p, iter = 5){
  Ahat = A
  for(i in 1:iter) {
    Ahat <- Ahat*Omega + (1- Omega) * svd_impute(Omega*Ahat/p, K)
  }
  Ahat
}

simple_impute <- function(A, Omega){
  diag(Omega) = 0
  phat = sum(A*Omega) / sum(Omega)
  diag(Omega) = 1
  Aimp = A*Omega + phat*(1-Omega)
  Aimp
}

loss_L1 <- function(A, Omega, Ahat) {
  diag(Omega) = 1
  Ahat = pmin(pmax(0, Ahat),1)
  sum( abs((A-Ahat)*(1-Omega))/sum(1-Omega))
}

loss_L2 <- function(A, Omega, Ahat) {
  diag(Omega) = 1
  Ahat = pmin(pmax(0, Ahat),1)
  sum( abs((A-Ahat)*(1-Omega))^2 /sum(1-Omega))
}

loss_binomialdeviance <- function(A, Omega, Ahat) {
  diag(Omega) = 1
  fold <- which(Omega == 0 & upper.tri(Omega), arr.ind = TRUE)
  A_w <- A[fold]
  Ahat_w <- Ahat[fold]
  Ahat_w <- pmin(pmax(Ahat_w, 1e-2), 1-1e-2)
  -sum(A_w * log(Ahat_w) +(1-A_w) *log(1-Ahat_w)) / length(Ahat_w)
}


sample_nodepairs <- function(N, p) {
  selected = matrix(0, N, N)
  selected[upper.tri(selected)] = rbinom(N*(N-1)/2, size = 1, p)
  selected <- selected + t(selected)
  selected
}

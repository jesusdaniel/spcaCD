library("MASS")

hard_thresh = function(A, g_threshold){
  A[(A) <= g_threshold] = 0
  return(A)
}
soft_thresh = function(A, g_threshold){
  return(sign(A)*pmax((A)-g_threshold,0))
}


#SPCA by thresholding (Ma, 2013)
spca_qr <- function(A, K, gamma = 0, lambdavec = rep(1, K), Z0 = NULL, thresh = soft_thresh,
                    num_iter = 50, TOL = 1e-04, save_path = T) {
  if(is.null(lambdavec))
    lambdavec <- abs(eigen(A)$values[1:K])
  lambda_mat =  matrix(kronecker(gamma * lambdavec, rep(1, nrow(A))), nrow = nrow(A))
  if(is.null(Z0))
    Z0 = initialization(ncol(A), K)
  Z0 = qr.Q(qr(Z0))
  Z = Z0                  #initial value
  tolerance = Inf;     iter = 1           #convergence
  tol_path = c();      Z_path = list()    #save path
  Zest = Z
  Z_path[[iter]] = Z0
  while(iter <= num_iter & tolerance >= TOL) {
    U = A %*%Z
    Y = thresh(as.matrix(U), 0*lambda_mat)
    QR = qr(Y)
    Z <- qr.Q(QR)

    iter = iter + 1
    tolerance = tolerance_function(Zest,Z)
    tol_path = c(tol_path, tolerance)
    if(save_path)
      Z_path[[iter]] = Z
    Zest = Z
  }
  return(list(Z = as.matrix(Z), iterations = iter, tolerance = tolerance,
              tol_path = tol_path, Z_path = Z_path))
}


spca_Z <- function(A, K, lambda = 0.9, lnorm = 1,
                   Z0 = NULL, num_iter = 50, TOL = 1e-04, save_path = T,
                   thresh = hard_thresh) {
  #browser()
  lambda_mat = array(lambda, dim = c(nrow(A), K))#threshold matrix
  lnorm_exp = 1/lnorm;
  if(is.null(Z0))
    Z0 = initialization(ncol(A), K)
  Z0 = t(apply(Z0, 1,function(x) x/(sum(abs(x)^lnorm)+1*(sum(abs(x))<1e-16))^lnorm_exp))

  Z = Z0                  #initial value
  tolerance = Inf;     iter = 1           #convergence
  tol_path = c();      Z_path = list()    #save path

  Zest = Z
  Z_path[[iter]] = Z0
  while(iter <= num_iter & tolerance >= TOL) {
    # Multiplication
    U = A %*%Z
    #Column normalization
    V = U%*%diag(apply(Z, 2,function(x) 1/(sum(abs(x)^lnorm)+1*(sum(abs(x))<1e-16))^lnorm_exp))
    #V = U %*% ginv(as.matrix(crossprod(Z)))
    # Row normalization by sup norm for thresholding
    V = t(apply(V, 1,function(x) x/(max(abs(x))+1*(sum(abs(x))<1e-16))))
    #Thresholding
    Y = thresh(as.matrix(V), lambda_mat)
    #Row renormalization - lnorm
    Z = t(apply(Y, 1,function(x) x/(sum(abs(x)^lnorm)+1*(sum(abs(x))<1e-16))^lnorm_exp))
    iter = iter + 1
    tolerance = tolerance_function(Zest,Z)
    tol_path = c(tol_path, tolerance)
    if(save_path)
      Z_path[[iter]] = Z
    Zest = Z
  }
  return(list(Z = as.matrix(Z), iterations = iter, tolerance = tolerance,
              tol_path = tol_path, Z_path = Z_path))
}



spca_DCZ <- function(A, K, lambda = 0.9, lnorm = 2,
                     Z0 = NULL, num_iter = 50, TOL = 1e-08, save_path = F,
                     thresh = hard_thresh) {
  N = nrow(A)
  lnorm_exp = 1/lnorm;
  if(is.null(Z0))
    Z0 = initialization(ncol(A), K)
  #Z0 = apply(Z0, 2,function(x) x/(sum(abs(x)^lnorm)+1*(sum(abs(x))<1e-16))^lnorm_exp)
  Z0 = apply(Z0, 2,function(x) x/(sum(abs(x))+1*(sum(abs(x))<1e-16)))
  Z = Z0                  #initial value
  tolerance = Inf;     iter = 1           #convergence
  tol_path = c();      Z_path = list()    #save path

  Zest = Z
  Z_path[[iter]] = Z0
  while(iter <= num_iter & tolerance >= TOL) {
    U = A %*%Z
    #V = U%*%diag(apply(Z, 2,function(x) 1/(sum(abs(x)^lnorm)+1*(sum(abs(x))<1e-16))^lnorm_exp))
    #V = U %*% ginv(as.matrix(crossprod(Z)))
    Lambdatemp = matrix(lambda*apply(U, 1,function(x) (max(abs(x))+1*(sum(abs(x))<1e-16))), ncol = K, nrow = N)
    Y = thresh(as.matrix(U), Lambdatemp)
    #Z = apply(Y, 2,function(x) x/(sum(abs(x)^lnorm)+1*(sum(abs(x))<1e-16))^lnorm_exp)
    Z = apply(Y, 2,function(x) x/(sum(abs(x))+1*(sum(abs(x))<1e-16)))
    iter = iter + 1
    tolerance = tolerance_function(Zest,Z)
    tol_path = c(tol_path, tolerance)
    if(save_path)
      Z_path[[iter]] = Z
    Zest = Z
  }
  return(list(Z = as.matrix(Z), iterations = iter, tolerance = tolerance,
              tol_path = tol_path, Z_path = Z_path))
}

spca_DCZ_inv <- function(A, K, lambda = 0.9, lnorm = 2,
                     Z0 = NULL, num_iter = 50, TOL = 1e-08, save_path = F,
                     thresh = hard_thresh) {
  N = nrow(A)
  lnorm_exp = 1/lnorm;
  if(is.null(Z0))
    Z0 = initialization(ncol(A), K)
  Z0 = apply(Z0, 2,function(x) x/(sum(abs(x)^lnorm)+1*(sum(abs(x))<1e-16))^lnorm_exp)

  Z = Z0                  #initial value
  tolerance = Inf;     iter = 1           #convergence
  tol_path = c();      Z_path = list()    #save path

  Zest = Z
  Z_path[[iter]] = Z0
  while(iter <= num_iter & tolerance >= TOL) {
    U = A %*%Z
    #V = U%*%diag(apply(Z, 2,function(x) 1/(sum(abs(x)^lnorm)+1*(sum(abs(x))<1e-16))^lnorm_exp))
    V = U %*% ginv(as.matrix(crossprod(Z)))
    Lambdatemp = matrix(lambda*apply(V, 1,function(x) (max(abs(x))+1*(sum(abs(x))<1e-16))), ncol = K, nrow = N)
    Y = thresh(as.matrix(V), Lambdatemp)
    Z = apply(Y, 2,function(x) x/(sum(abs(x)^lnorm)+1*(sum(abs(x))<1e-16))^lnorm_exp)
    iter = iter + 1
    tolerance = tolerance_function(Zest,Z)
    tol_path = c(tol_path, tolerance)
    if(save_path)
      Z_path[[iter]] = Z
    Zest = Z
  }
  return(list(Z = as.matrix(Z), iterations = iter, tolerance = tolerance,
              tol_path = tol_path, Z_path = Z_path))
}


spca_eigenbasis <- function(A, K, lambda = 0.3, lnorm = 2,
                              Z0 = NULL, num_iter=30, TOL = 1e-04, save_path = T,
                            thresh = hard_thresh) {
  N = nrow(A)
  lnorm_exp = 1/lnorm;
  if(is.null(Z0))
    Z0 = initialization(ncol(A), K)
  Z0 = apply(Z0, 2,function(x) x/(sum(abs(x)^lnorm)+1*(sum(abs(x))<1e-16))^lnorm_exp)

  Z = Z0                  #initial value
  tolerance = Inf;     iter = 1           #convergence
  tol_path = c();      Z_path = list()    #save path

  Zest = Z
  Z_path[[iter]] = Z0

  while(iter <= num_iter & tolerance >= TOL) {
    AZ = A %*%Z
    ZAZ = crossprod(Z, AZ)
    ZZ = crossprod(Z)
    V = AZ %*% (ginv(as.matrix(ZAZ)) %*% ZZ)
    Lambdatemp = matrix(lambda * apply(abs(V), 1,function(x) max(x)), nrow = N, ncol = K)
    Y = thresh(as.matrix(V), Lambdatemp)
    #Z = apply(Y, 2,function(x) x/(sum(abs(x)^lnorm)+1*(sum(abs(x))<1e-16))^lnorm_exp)
    Z = apply(Y, 2,function(x) x/(sum(abs(x)^2)+1*(sum(abs(x))<1e-16))^0.5)
    iter = iter + 1
    tolerance = tolerance_function(Zest,Z)
    tol_path = c(tol_path, tolerance)
    if(save_path)
      Z_path[[iter]] = Z
    Zest = Z
  }
  return(list(Z = as.matrix(Z), iterations = iter, tolerance = tolerance,
              tol_path = tol_path, Z_path = Z_path))
}


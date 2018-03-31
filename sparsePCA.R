library("MASS")

hard_thresh = function(A, g_threshold){
  A[abs(A) <= g_threshold] = 0
  return(A)
}
soft_thresh = function(A, g_threshold){
  return(sign(A)*pmax(abs(A)-g_threshold,0))
}



spca_iterative <- function(A, K, lambda, lnorm = 1, 
                           Z0, num_iter, TOL = 1e-04, save_path = T,
                           thresh = soft_thresh) {
  #browser()
  lambda_mat = array(lambda, dim = c(nrow(A), K))#threshold matrix
  lnorm_exp = 1/lnorm;
  Z0 = t(apply(Z0, 1,function(x) x/(sum(abs(x)^lnorm)+1*(sum(abs(x))<1e-16))^lnorm_exp))
  Z = Z0                  #initial value
  tolerance = Inf;     iter = 1           #convergence
  tol_path = c();      Z_path = list()    #save path
  
  Zest = Z
  Z_path[[iter]] = Z0
  while(iter <= num_iter & tolerance >= TOL) {
    
    U = A %*%Z
    #Normalize column
    #Z = U%*%diag(apply(Z, 2,function(x) 1/(sum(abs(x)^lnorm)+1*(sum(abs(x))<1e-16))^lnorm_exp))
    Z = U%*%solve(crossprod(Z))
    # Row normalization
    Z = t(apply(Z, 1,function(x) x/(max(abs(x))+1*(sum(abs(x))<1e-16))))
    Z = thresh(as.matrix(Z), lambda_mat)
    Z = t(apply(Z, 1,function(x) x/(sum(abs(x)^lnorm)+1*(sum(abs(x))<1e-16))^lnorm_exp))
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


spca_iterative_DC <- function(A, K, lambda, lnorm, 
                           Z0, num_iter, TOL = 1e-08, save_path = T) {
  
  #browser()
  lambda_mat = array(lambda, dim = c(nrow(A), K))#threshold matrix
  lnorm_exp = 1/lnorm;
  Z0 = t(apply(Z0, 1,function(x) x/(sum(abs(x)^lnorm)+1*(sum(abs(x))<1e-16))^lnorm_exp))
  Z = Z0                  #initial value
  tolerance = Inf;     iter = 1           #convergence
  tol_path = c();      Z_path = list()    #save path
  
  Zest = Z
  Z_path[[iter]] = Z0
  while(iter <= num_iter & tolerance >= TOL) {
    #Normalize column
    U = A %*%Z
    #Z = U%*%diag(apply(Z, 2,function(x) 1/(sum(abs(x)^lnorm)+1*(sum(abs(x))<1e-16))^lnorm_exp))
    Z = U%*% solve(crossprod(Z))
    #Z = t(apply(Z, 1,function(x) x/(max(abs(x))+1*(sum(abs(x))<1e-16))))
    theta_hat = apply(Z,1,max)
    Z = soft_thresh(as.matrix(Z), theta_hat*lambda_mat)
    theta_hat = apply(Z,1,sum)
    theta_hat = theta_hat/sum(theta_hat)
    Z = t(apply(Z, 1,function(x) x/(sum(abs(x)^lnorm)+1*(sum(abs(x))<1e-16))^lnorm_exp))
    Z = theta_hat*Z
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

tolerance_function = function(u,v) {
  #nonzeros = sum(u -  v %*% diag(sign(diag(crossprod(u,v))))!=0)
  #if(nonzeros==0)
  #  return(0)
  #norm(u -  v %*% diag(sign(diag(crossprod(u,v)))))/sqrt(nonzeros)
  norm(u-v)/sum(u!=0)
}




spca_iterative_DC2 <- function(A, K, lambda, lnorm, thresh = hard_thresh,
                           Z0, num_iter, TOL = 1e-04, save_path = T) {
  #browser()
  lambda_mat = array(lambda, dim = c(nrow(A), K))#threshold matrix
  lnorm_exp = 1/lnorm;
  Z0 = apply(Z0, 2,function(x) x/(sum(abs(x)^lnorm)+1*(sum(abs(x))<1e-16))^lnorm_exp)
  Z = Z0                  #initial value
  tolerance = Inf;     iter = 1           #convergence
  tol_path = c();      Z_path = list()    #save path
  
  Zest = Z
  Z_path[[iter]] = Z0
  while(iter <= num_iter & tolerance >= TOL) {
    #Normalize column
    
    projeig <- eig %*% solve(crossprod(eig)) %*% t(eig)
    v1 <- abs(eigen(A[1:(N/K), 1:(N/K)])$vectors[,1])
    v2 <- abs(eigen(A[(N/K+1):N,(N/K+1):N])$vectors[,1])
    Vs = cbind(c(v1,rep(0, N/K)), c(rep(0,N/K), v2))
    
    #errorsLists
    error_thresh_mult <- list()
    error_mult_eig <- list()
    error_Z_eig <- list()
    error_Z_Ztar <- list()
    error_Z0_Zstar <- list()
    error_V_Zstar <- list()
    zeros <- list()
    i = 1
    
    
    projZ0 = Z %*% solve(crossprod(Z)) %*% t(Z)
    error_Z0_Zstar[[i]] = norm(projZ0 - projZstar, "2")
    U = A %*%Z
    #table(U[,1]>U[,2],Ztrue)
    #sapply(1:2, function(i) sapply(1:2, function(j) 
    #  norm(tcrossprod(Vs[,i])/as.double(crossprod(Vs[,i])) -  tcrossprod(U[,j])/as.double(crossprod(U[,j])), "2")))
    V = U%*%diag(apply(U, 2,function(x) 1/(sum(abs(x)^lnorm)+1*(sum(abs(x))<1e-16))^lnorm_exp))
    #plot(V[,1])
    #points(V[,2], col = "red")
    
    Lambdatemp = diag(apply(V, 1,function(x) (max(abs(x))+1*(sum(abs(x))<1e-16)))) %*% lambda_mat
    Y = soft_thresh(as.matrix(V), Lambdatemp)
    par(mfrow=c(1,2))
    plot(V[,1], V[,2], col = Ztrue, xlim = c(0, max(V[,1])), ylim = c(0,max(V[,2])))
    plot(Y[,1], Y[,2], col = Ztrue, xlim = c(0, max(Y[,1])), ylim = c(0,max(Y[,2])))
    Z = Y %*%diag(apply(Y, 2,function(x) 1/(sum(abs(x)^lnorm)+1*(sum(abs(x))<1e-16))^lnorm_exp))
    #sapply(1:2, function(i) sapply(1:2, function(j) 
    #  norm(tcrossprod(Vs[,i])/as.double(crossprod(Vs[,i])) -  tcrossprod(Z[,j])/as.double(crossprod(Z[,j])), "2")))
    #sapply(1:2, function(i) sapply(1:2, function(j) 
    #  sqrt(sum((Vs[,i] - Z[,j])^2))))
    #Z
    projZ = Z %*% solve(crossprod(Z)) %*% t(Z)
    projV = V %*% solve(crossprod(V)) %*% t(V)
    error_thresh_mult[[i]] = norm(projV - projZ, "2")
    #error_mult_eig[[i]] = norm(projV - projeig, "2")
    #error_Z_eig[[i]] = norm(projZ - projeig, "2")
    error_Z_Ztar[[i]] = norm(projZ - projZstar, "2")
    error_V_Zstar[[i]] = norm(projV - projZstar, "2")
    #zeros[[i]] = sum(Z==0)
    i = i+1
    
    
    plot(unlist(error_Z0_Zstar))
    points(unlist(error_Z_Ztar), col = "red")
    points(unlist(error_V_Zstar) + unlist(error_thresh_mult),
          col = "blue")
    
    plot(unlist(error_mult_eig))
    plot(unlist(error_Z_eig))
    plot(unlist(error_Z_Ztar))
    plot(unlist(zeros))
    
    Zstar = Z
    projZstar = Zstar %*% solve(crossprod(Zstar)) %*% t(Zstar)
    norm(projeig - projZstar, "2")
    
    
    plot(eig[,1], eig[,2], col = Ztrue, xlim = c(min(eig[,1]),0))
    plot(V[,1], V[,2], col = Ztrue, xlim = c(0, max(V[,1])), ylim = c(0,max(V[,2])))
    
    u1ind = which(Zstar[,1]>0)
    Asub1 <- A[u1ind, u1ind]
    u <- abs(eigen(Asub1)$vectors[,1])
    plot(v1, Zstar[u1ind,1])
    u2ind = which(Zstar[,2]>0)
    Asub2 <- A[u2ind, u2ind]
    v2 <- abs(eigen(Asub2)$vectors[,1])
    plot(v2, Zstar[u2ind,2])
    #plot(v2- Zstar[u2ind,2])
    Z0u <- cbind(c(u, rep(0,100)), c(rep(0,100), v2))
    
    table(Ztrue, Z[,1])
    
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


spca_iterative_DC <- function(A, K, lambda, lnorm = 2, thresh = hard_thresh,
                              Z0, num_iter, TOL = 1e-04, save_path = T) {
  lambda_mat = array(lambda, dim = c(nrow(A), K))#threshold matrix
  lnorm_exp = 1/lnorm;
  Z0 = apply(Z0, 2,function(x) x/(sum(abs(x)^lnorm)+1*(sum(abs(x))<1e-16))^lnorm_exp)
  Z = Z0                  #initial value
  tolerance = Inf;     iter = 1           #convergence
  tol_path = c();      Z_path = list()    #save path
  
  Zest = Z
  Z_path[[iter]] = Z0
  while(iter <= num_iter & tolerance >= TOL) {
    U = A %*%Z
    #Normalize column
    #V = U%*%diag(apply(U, 2,function(x) 1/(sum(abs(x)^lnorm)+1*(sum(abs(x))<1e-16))^lnorm_exp))
    V = U%*%solve(crossprod(Z))
    Lambdatemp = diag(apply(V, 1,function(x) (max(abs(x))+1*(sum(abs(x))<1e-16)))) %*% lambda_mat
    Y = thresh(as.matrix(V), Lambdatemp)
    Z = Y %*%diag(apply(Y, 2,function(x) 1/(sum(abs(x)^lnorm)+1*(sum(abs(x))<1e-16))^lnorm_exp))
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


spca_iterative_inverse <- function(A, K, lambda, lnorm = 2, thresh = hard_thresh,
                              Z0, num_iter, TOL = 1e-04, save_path = T) {
  browser()
  lambda_mat = array(lambda, dim = c(nrow(A), K))#threshold matrix
  lnorm_exp = 1/lnorm;
  #Z0 = apply(Z0, 2,function(x) x/(sum(abs(x)^lnorm)+1*(sum(abs(x))<1e-16))^lnorm_exp)
  Z = Z0                 #initial value
  Zinvhat = solve(crossprod(Z))
  Bhat = Zinvhat %*% (crossprod(Z, A) %*% Z) %*% Zinvhat
  tolerance = Inf;     iter = 1           #convergence
  tol_path = c();      Z_path = list()    #save path
  
  Zest = Z
  Z_path[[iter]] = Z0
  while(iter <= num_iter & tolerance >= TOL) {
    U = A %*%Z
    #Normalize column
    Zinvhat = solve(crossprod(Z))
    Bhat = Zinvhat %*% (crossprod(U, Z)) %*% Zinvhat
    V = U%*%Zinvhat %*% solve(Bhat)
    Lambdatemp = diag(apply(V, 1,function(x) (max(abs(x))+1*(sum(abs(x))<1e-16)))) %*% lambda_mat
    Z = thresh(as.matrix(V), Lambdatemp)
    #Z = Z*(Z>0)
    Zinvhat = solve(crossprod(Z))
    Bhat = Zinvhat %*% (crossprod(U, Z)) %*% Zinvhat
    thetainv = sqrt(diag(Bhat))
    Z = as.matrix(Z %*% diag(thetainv))
    iter = iter + 1
    tolerance = tolerance_function(Zest,Z)
    tol_path = c(tol_path, tolerance)
    if(save_path)
      Z_path[[iter]] = Z
    Zest = Z
  }
  Zinvhat = solve(crossprod(Z))
  Bhat = Zinvhat %*% (crossprod(U, Z)) %*% Zinvhat
  thetainv = sqrt(diag(Bhat))
  return(list(Z = as.matrix(Z %*% diag(thetainv)), iterations = iter, tolerance = tolerance,
              tol_path = tol_path, Z_path = Z_path))
}


spca_iterative_other <- function(A, K, lambda, lnorm = 2, thresh = hard_thresh,
                                   Z0, num_iter, TOL = 1e-04, save_path = T) {
  browser()
  lambda_mat = array(lambda, dim = c(nrow(A), K))#threshold matrix
  lnorm_exp = 1/lnorm;
  #Z0 = apply(Z0, 2,function(x) x/(sum(abs(x)^lnorm)+1*(sum(abs(x))<1e-16))^lnorm_exp)
  Z = Z0                 #initial value
  Zinvhat = solve(crossprod(Z))
  Bhat = Zinvhat %*% (crossprod(Z, A) %*% Z) %*% Zinvhat
  tolerance = Inf;     iter = 1           #convergence
  tol_path = c();      Z_path = list()    #save path
  
  Zest = Z
  Z_path[[iter]] = Z0
  while(iter <= num_iter & tolerance >= TOL) {
    grad = -2*Z  %*% solve(crossprod(Z)) %*% (crossprod(Z,A) %*% Z) %*% solve(crossprod(Z))
    U = A %*%Z
    #Normalize column
    Zinvhat = solve(crossprod(Z))
    Bhat = Zinvhat %*% (crossprod(U, Z)) %*% Zinvhat
    V = U%*%Zinvhat %*% solve(Bhat)
    Lambdatemp = diag(apply(V, 1,function(x) (max(abs(x))+1*(sum(abs(x))<1e-16)))) %*% lambda_mat
    Z = thresh(as.matrix(V), Lambdatemp)
    #Z = Z*(Z>0)
    Zinvhat = solve(crossprod(Z))
    Bhat = Zinvhat %*% (crossprod(U, Z)) %*% Zinvhat
    thetainv = sqrt(diag(Bhat))
    Z = as.matrix(Z %*% diag(thetainv))
    iter = iter + 1
    tolerance = tolerance_function(Zest,Z)
    tol_path = c(tol_path, tolerance)
    if(save_path)
      Z_path[[iter]] = Z
    Zest = Z
  }
  Zinvhat = solve(crossprod(Z))
  Bhat = Zinvhat %*% (crossprod(U, Z)) %*% Zinvhat
  thetainv = sqrt(diag(Bhat))
  return(list(Z = as.matrix(Z %*% diag(thetainv)), iterations = iter, tolerance = tolerance,
              tol_path = tol_path, Z_path = Z_path))
}

spca_iterative_original <- function(A, K, lambda, lnorm, 
                               Z0, num_iter, TOL = 1e-04, save_path = T) {
  #browser()
  lambdas <- eigen(A)$values[1:K]
  lambda_mat =  matrix(kronecker(lambda * lambdas, rep(1, nrow(A))), nrow = nrow(A))
  lnorm_exp = 2;
  Z0 = apply(Z0, 2,function(x) x/(sum(abs(x)^lnorm)+1*(sum(abs(x))<1e-16))^(1/lnorm_exp))
  Z = Z0                  #initial value
  tolerance = Inf;     iter = 1           #convergence
  tol_path = c();      Z_path = list()    #save path
  
  Zest = Z
  Z_path[[iter]] = Z0
  while(iter <= num_iter & tolerance >= TOL) {
    #Normalize column
    
    U = A %*%Z
    #V = U%*%diag(1/sqrt(lambdas[1:2]))
    #diag(apply(Z, 2,function(x) 1/(sum(abs(x)^lnorm)+1*(sum(abs(x))<1e-16))^lnorm_exp))
    Y = soft_thresh(as.matrix(U), 0*lambda_mat)
    QR = qr(Y)
    Z <- qr.Q(QR)
    crossprod(Z,eig)
    #Y
    #Lambdatemp = diag(apply(V, 1,function(x) (max(abs(x))+1*(sum(abs(x))<1e-16)))) %*% lambda_mat
    #Y = soft_thresh(as.matrix(V), Lambdatemp)
    #Z = Y
    #Z
    projZ = Z %*% solve(crossprod(Z)) %*% t(Z)
    projV = V %*% solve(crossprod(V)) %*% t(V)
    norm(projZ - projV, "2")
    projeig <- eig %*% solve(crossprod(eig)) %*% t(eig)
    norm(projV - projeig, "2")
    
    plot(V[,1])
    plot(V[,2])
    Z=V
    
    plot(Z[,1])
    plot(Z[,2])
    table(Ztrue, Z[,1])
    
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




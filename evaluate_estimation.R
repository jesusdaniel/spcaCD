# Measure accuracy
table2x2 <- function(G1, G2) {
  #browser()
  tab = table(G1,G2)
  if(length(tab) == 4){
    return(tab)
  }
  tab2 = array(0, dim = c(2,2))
  if(length(tab) == 1){
    tab2[unique(G1)+1,unique(G2)+1] = tab
    return(tab2)
  }
  if(length(unique(G1)) == 1) {
    tab2[unique(G1)+1,] = tab
  }else{
    tab2[,unique(G2)+1] = tab
  }
  return(tab2)
}
binarize <- function(Zest) {
  t(apply(Zest,1,function(x){u = 1*(abs(x)==max(abs(x))); if(sum(u)>1) return(0*u); return(u)}))
}
accuracy <- function(Zest, Z) {
  K = ncol(Z)
  maxaccur = 0
  for(permut in combinat::permn(K)) {
    accur = sum(Zest[,permut]==Z &Z==1)/nrow(Z)
    if(accur> maxaccur)
      maxaccur = accur
  }
  return(maxaccur)
}
accuracy_osbm_2k <- function(Zest, Z) {
  K = 2
  if(sum(diag(crossprod(Zest,Z)))>sum(diag(crossprod(Zest[,c(2,1)],Z)))) {
    communities = c(1,2)
  }else{communities = c(2,1)}
  tot = 0
  Zest = 1*(Zest!=0)
  Z = 1*(Z!=0)
  for(i in 1:K)
    tot = tot + sum(Zest[,which(communities==i)]==Z[,i])
  return(tot/(2*nrow(Z)))
}
# L2 error
l2error <- function(Zt,Zest) {
  return(min(sapply(combinat::permn(1:K), function(sigm) {
    sum(sapply(1:K, function(j) sum((Zt[,sigm[j]] - Zest[,j])^2)))
    })))
}

trueposfalsepos_continuous <- function(Zest, Z) {
  Zbin <- binarize(Zest)
  Ztruebin = 1*(Z!=0)
  best = which.max(sapply(permn(1:ncol(Zbin)), function(u) sum(diag(crossprod(Zbin[,u], Ztruebin)))))
  Zest = Zest[, permn(1:ncol(Zest))[[best]]]
  truepos = c()
  falsepos = c()
  sparsity = c()
  for(lambda in seq(0.05,1, 0.01)){
    Zbin = 1*(abs(Zest)>=lambda)
    truepos = c(truepos, sum(Zbin*Ztruebin))
    sparsity = c(sparsity, sum(Zbin!=0))
    falsepos = c(falsepos, sum(Zbin*(1-Ztruebin)))
  }
  plot(falsepos, truepos)
}
trueposfalsepos_path <- function(Zpath, Z) {
  Ztruebin = 1*(Z!=0)
  truepos = c()
  falsepos = c()
  sparsity = c()
  for(Zest in Zpath) {
    Zbin <- binarize(Zest)
    best = which.max(sapply(permn(1:ncol(Zbin)), function(u) sum(diag(crossprod(Zbin[,u], Ztruebin)))))
    Zest = Zest[, permn(1:ncol(Zest))[[best]]]
    Zbin = 1*(abs(Zest)>0)
    truepos = c(truepos, sum(Zbin*Ztruebin))
    sparsity = c(sparsity, sum(Zbin!=0))
    falsepos = c(falsepos, sum(Zbin*(1-Ztruebin)))
  }
  points(falsepos, truepos, col = "red")
}
Zest = sol_occam
Zest = sol_mSCORE
exniv_path <- function(Zpath, Z) {
  exniv = c()
  sparsity = c()
  for(Zest in Zpath){
    exniv = c(exniv, exNiV((Zest!=0)*1, Z))
    sparsity = c(sparsity, sum(Zest!=0))
  }
  points(sparsity, exniv, col = "red")
}
exniv_sparsity <- function(Zest, Z) {
  exniv = c()
  sparsity = c()
  for(lambda in seq(0.05,1, 0.01)){
    exniv = c(exniv, exNiV((Zest>lambda), Z))
    sparsity = c(sparsity, sum(Zest>lambda))
  }
  plot(sparsity, exniv, col = "blue")
}
best_exNIV_threshold <- function(Zest, Z) {
  Z = Z>0
  bestEXNIV = 0
  for(lambda in seq(0,0.99,0.01))
    bestEXNIV = max(bestEXNIV, (exNiV((Zest>lambda), Z)))
  return(bestEXNIV)
}
exNIV_threshold <- function(Zest, Z, Thr = NULL) {
  if(is.null(Thr))
    Thr = 1/sqrt(ncol(Z))
  exNiV((Zest>Thr), Z>0)
}
exNiV <- function(Zest, Z) {
  Z = abs(sign(Z))
  n = nrow(Z); K = ncol(Z)
  P_hat = apply(Zest,2,sum)/n;   P = apply(Z,2,sum)/n
  H_hat = -(P_hat*log(P_hat) +(1-P_hat)*log(1-P_hat))
  H_hat = ifelse(is.nan(H_hat),0,H_hat)
  H = -(P*log(P) +(1-P)*log(1-P))
  H = ifelse(is.nan(H),0,H)
  
  H_Z_Zest = -apply(Zest,2, function(x) apply(Z,2,table2x2,x))/n*log(apply(Zest,2, function(x) apply(Z,2,table2x2,x))/n)
  H_Z_Zest = ifelse(is.nan(H_Z_Zest),0,H_Z_Zest)
  H_Z_Zest = t(kronecker(diag(K),rep(1,4))) %*% H_Z_Zest
  
  H_Zest_Z = -apply(Z,2, function(x) apply(Zest,2,table2x2,x))/n*log(apply(Z,2, function(x) apply(Zest,2,table2x2,x))/n)
  H_Zest_Z = ifelse(is.nan(H_Zest_Z),0,H_Zest_Z)
  H_Zest_Z  = t(kronecker(diag(K),rep(1,4))) %*% H_Zest_Z
  
  H_Zest_giv_Z = H_Z_Zest - H
  Hbar_Zest_giv_Z = H_Zest_giv_Z/H_hat
  Hbar_Zest_giv_Z = ifelse(is.nan(Hbar_Zest_giv_Z)|Hbar_Zest_giv_Z==-Inf,0,Hbar_Zest_giv_Z)
  
  H_Z_giv_Zest = H_Zest_Z - H_hat
  Hbar_Z_giv_Zest = H_Z_giv_Zest/H
  Hbar_Z_giv_Zest = ifelse(is.nan(Hbar_Z_giv_Zest),0,Hbar_Z_giv_Zest)
  
  return(1-min(sapply(combinat::permn(1:K), function(sigm){
    sum(sapply(1:K, function(j) {
      Hbar_Zest_giv_Z[sigm[j],j] + Hbar_Z_giv_Zest[j,sigm[j]]
    }))
  }))/(2*K))
}

loglikelihood <- function(A,Z,B) {
  W = tcrossprod(tcrossprod(Z, B),Z)
  L = A*log(W) + (1-A)*log(1-W)
  return(sum(L[upper.tri(L,diag = T)]))
}



find_purenodes <- function(Zest,com,K) {
  max_com = K
  which(apply(Zest[,(1:max_com)[-com]],1,sum)==0 & Zest[,com]!=0)
}

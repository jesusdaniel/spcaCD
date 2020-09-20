Z_from_communities <- function(trueC) {
  K = length(unique(trueC))
  Z = matrix(0, length(trueC), K)
  for(i in 1:K) {
    Z[trueC==i, i] = 1
  }
  return(Z)
}
# Calculate Normalized variation of information
# Z, Zest
exNVI <- function(Zest, Z) {
  Z = abs(sign(Z))
  Zest = abs(sign(Zest))
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



# For continuous memberships
exNVI_threshold <- function(Zest, Z, Thr = NULL) {
  if(is.null(Thr))
    Thr = 1/sqrt(ncol(Z))
  exNVI((Zest>Thr), Z>0)
}

# Best threshold for continuous memberships
best_exNVI_threshold <- function(Zest, Z) {
  Z = Z>0
  bestEXNIV = 0
  for(lambda in seq(0,0.99,0.03))
    bestEXNIV = max(bestEXNIV, (exNVI((Zest>lambda), Z)))
  return(bestEXNIV)
}

exNVI_continuouspath <- function(Zest, Z, lambdapath =  seq(0,0.99,0.03)) {
  Z = Z>0
  acc <- sapply(lambdapath, function(lambda) exNVI((Zest>lambda), Z))
  sp <- nrow(Z)*ncol(Z) - sapply(lambdapath, function(lambda) sum(Zest>lambda))
  return(list(acc = acc, sp = sp))
}

# Distance between subspaces
distance_subs <- function(Z1, Z2) {
  if(ncol(Z1) != ncol(Z2))
    1
  base::norm(as.matrix(Z1 %*% ginv(as.matrix(crossprod(Z1))) %*% t(Z1) - 
                         Z2 %*% ginv(as.matrix(crossprod(Z2))) %*% t(Z2)), "2")
}


distance_L2 <- function(A, Z) {
  norm(A  - Z %*% (ginv(as.matrix(crossprod(Z))) %*% 
         ((crossprod(Z,A) %*% Z) %*% ginv(as.matrix(crossprod(Z))))) %*% t(Z), "F")
}

# non-ovelapping memberships accuracy
binarize <- function(Zest) {
  t(apply(Zest,1,function(x){u = 1*(abs(x)==max(abs(x))); if(sum(u)>1) return(0*u); return(u)}))
}

binarize_one <- function(Zest) {
  apply(abs(Zest), 1, which.max)
}

accuracy <- function(Zest, Z) {
  nonoverlap = which(apply(Z!=0,1, sum) ==1)
  Zest = binarize(Zest[nonoverlap,])
  Z = binarize(Z[nonoverlap,])
  K <- ncol(Z)
  maxaccur = 0
  for(permut in combinat::permn(K)) {
    accur = sum(Zest[,permut]==Z &Z==1)/nrow(Z)
    if(accur> maxaccur)
      maxaccur = accur
  }
  return(maxaccur)
}



# overlapping accuracy -------------------------------------
overlapping_accuracy <- function(Zest, Z) {
  Zest = 1*(Zest !=0 )
  Z = 1*(Z != 0)
  K <- ncol(Z)
  maxaccur = 0
  for(permut in combinat::permn(K)) {
    accur = sum(Zest[,permut]==Z)/(nrow(Z) * K)
    if(accur> maxaccur)
      maxaccur = accur
  }
  return(maxaccur)
}
best_overlapaccuracy <- function(Zest, Z) {
  Z = Z>0
  bestoa = 0
  for(lambda in seq(0,0.99,0.01))
    bestoa = max(bestoa, (overlapping_accuracy((Zest>lambda), Z)))
  return(bestoa)
}


trueposfalseposrate_continuous <- function(Zest, Z) {
  Zbin <- binarize(Zest)
  Ztruebin = 1*(Z!=0)
  best = which.max(sapply(permn(1:ncol(Zbin)), function(u) sum(diag(crossprod(Zbin[,u], Ztruebin)))))
  Zest = Zest[, permn(1:ncol(Zest))[[best]]]
  Zest = t(apply(Zest, 1, function(x) x/(max(x)+1*(max(x)==0))))
  truepos = c()
  falsepos = c()
  sparsity = c()
  for(lambda in seq(0.05,1, 0.01)){
    Zbin = 1*(abs(Zest)>=lambda)
    truepos = c(truepos, sum(Zbin*Ztruebin))
    sparsity = c(sparsity, sum(Zbin!=0))
    falsepos = c(falsepos, sum(Zbin*(1-Ztruebin)))
  }
  total_truepos = sum(Z!=0)
  total_falsepos = sum(Z==0)
  return(list(falsepos = falsepos/total_falsepos, truepos = truepos/total_truepos))
}

trueposfalseposrate_path <- function(Zpath, Z) {
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
  total_truepos = sum(Z!=0)
  total_falsepos = sum(Z==0)
  return(list(falsepos = falsepos/total_falsepos, truepos = truepos/total_truepos))
}

trueposfalseposrate_overlap_continuous <- function(Zest, Z) {
  Zbin <- binarize(Zest)
  Ztruebin = 1*(Z!=0)
  Overlap_true <- 1*(apply(Ztruebin, 1, sum) > 1)
  best = which.max(sapply(permn(1:ncol(Zbin)), function(u) sum(diag(crossprod(Zbin[,u], Ztruebin)))))
  Zest = Zest[, permn(1:ncol(Zest))[[best]]]
  Zest = t(apply(Zest, 1, function(x) x/(max(x)+1*(max(x)==0))))
  truepos = c()
  falsepos = c()
  for(lambda in seq(0.05,1, 0.01)){
    Zbin = 1*(abs(Zest)>=lambda)
    Overlap_hat = 1*(apply(Zbin, 1, sum) > 1)
    truepos = c(truepos, sum(Overlap_hat*Overlap_true))
    falsepos = c(falsepos, sum(Overlap_hat*(1-Overlap_true)))
  }
  total_truepos = sum(Overlap_true!=0)
  total_falsepos = sum(Overlap_true==0)
  return(list(falsepos = falsepos/total_falsepos, truepos = truepos/total_truepos))
}

trueposfalseposrate_overlap_path <- function(Zpath, Z) {
  Ztruebin = 1*(Z!=0)
  Overlap_true <- 1*(apply(Ztruebin, 1, sum) > 1)
  truepos = c()
  falsepos = c()
  sparsity = c()
  for(Zest in Zpath) {
    Zbin <- binarize(Zest)
    best = which.max(sapply(permn(1:ncol(Zbin)), function(u) sum(diag(crossprod(Zbin[,u], Ztruebin)))))
    Zest = Zest[, permn(1:ncol(Zest))[[best]]]
    Zbin = 1*(abs(Zest)>0)
    Overlap_hat = 1*(apply(Zbin, 1, sum) > 1)
    truepos = c(truepos, sum(Overlap_hat*Overlap_true))
    falsepos = c(falsepos, sum(Overlap_hat*(1-Overlap_true)))
  }
  total_truepos = sum(Overlap_true!=0)
  total_falsepos = sum(Overlap_true==0)
  return(list(falsepos = falsepos/total_falsepos, truepos = truepos/total_truepos))
}


trapezoid <- function(x, y) {
  # numerical integral of fun from a to b
  # using the trapezoid rule with n subdivisions
  # assume a < b and n is a positive integer
  h <- diff(x)
  n <- length(y) - 1
  s <- sum(y[1:n] *h + (y[2:(n+1)] - y[1:n]) *h/2)
  return(s)
}
AUC_ROC <- function(Zpath, Z) {
  tfpos <- trueposfalseposrate_path(Zpath, Z)
  x = c(0, sort(tfpos$falsepos), 1)
  y = c(0, tfpos$truepos[order(tfpos$falsepos)], 1)
  return(trapezoid(x, y ))
}

AUC_ROC_continuous <- function(Zest, Z) {
  tfpos <- trueposfalseposrate_continuous(Zest, Z)
  x = c(0, sort(tfpos$falsepos), 1)
  y = c(0, tfpos$truepos[order(tfpos$falsepos)], 1)
  return(trapezoid(x, y ))
}

AUC_ROC_overlap <- function(Zpath, Z) {
  tfpos <- trueposfalseposrate_overlap_path(Zpath, Z)
  x = c(0, sort(tfpos$falsepos), 1)
  y = c(0, tfpos$truepos[order(tfpos$falsepos)], 1)
  return(trapezoid(x, y ))
}

AUC_ROC_overlap_continuous <- function(Zest, Z) {
  tfpos <- trueposfalseposrate_overlap_continuous(Zest, Z)
  x = c(0, sort(tfpos$falsepos), 1)
  y = c(0, tfpos$truepos[order(tfpos$falsepos)], 1)
  return(trapezoid(x, y ))
}

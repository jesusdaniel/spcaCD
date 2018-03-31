library(gtools)
## Functions
## Geometric functions
cut_f = function(x, thr){
  x[x<thr]=0
  return(x/sum(x))
}
barc = function(beta, doc, w=1){
  beta = rbind(t(beta), rep(1, K))
  f = lsfit(beta, c(doc,1), wt=c(rep(1,V),w), intercept=FALSE)
  return(c(f$coef))
}

convert = function(beta){
  if (nrow(beta)==2) return(beta[2,] - beta[1,])
  return(apply(beta[-1,], 1, function(x) x - beta[1,]))
}

dist_to_s = function(doc, beta){
  if (is.vector(beta)) return(list(d=sum((beta - doc) ^ 2), f=TRUE))
  coord = lsfit(convert(beta), doc-beta[1,], intercept=FALSE)
  if (all(coord$coef>=0) && sum(coord$coef)<=1) return(list(d = sum(coord$residuals^2), f=TRUE)) # numeric issues with 1?
  if (any(coord$coef<0)) return(list(ind = c(1,which(coord$coef>0)+1), f=FALSE))
  if (all(coord$coef>=0) && sum(coord$coef)>1) return(list(ind=2:nrow(beta), f=FALSE))
  print('error')
}

get_dist = function(doc, beta){
  f = FALSE
  repeat{
    out = dist_to_s(doc, beta)
    if (out$f) break
    beta = beta[out$ind,]
  }
  return(out$d)
}

proj_on_s = function(doc, beta, ind){
  beta = beta[ind,]
  crd = rep(0,K)
  if (is.vector(beta)){
    crd[ind] = 1
    return(list(c=crd, f=TRUE))
    print('caution')
  }
  coord = lsfit(convert(beta), doc-beta[1,], intercept=FALSE)
  if (all(coord$coef>=0) && sum(coord$coef)<=1){
    crd[ind] = c(1-sum(coord$coef), coord$coef)
    return(list(c = crd, f=TRUE))
  }
  if (any(coord$coef<0)) return(list(ind = ind[c(1,which(coord$coef>0)+1)], f=FALSE))
  if (all(coord$coef>=0) && sum(coord$coef)>1) return(list(ind=ind[-1], f=FALSE))
  print('error')
}

get_coord = function(doc, beta){
  f = FALSE
  ind = 1:K
  repeat{
    out = proj_on_s(doc, beta, ind)
    if (out$f) break
    ind = out$ind
  }
  return(out$c)
}

euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))

simul = function(V, K, N, L, eta, alpha, H){
  beta_t = rdirichlet(K, rep(eta, V)) #no cut
  teta_t = rdirichlet(N+H, rep(alpha, K))
  simplex = teta_t %*% beta_t
  ## generate words
  wdf = t(apply(simplex, 1, function(x) rmultinom(1, L, x))) 
  return(list(wdf=wdf[1:N,], beta_t=beta_t, test=wdf[(N+1):(N+H),]))
}

## Fixed
# N = 3000
# V = 1200
# L = 1000
# H = 100
# K = 5
# eta = 0.1
# alpha = 0.1
# 
# s = simul(V, K, N, L, eta, alpha, H)
# wdf = s$wdf
# wdfn = t(apply(wdf, 1, function(x) x/sum(x)))
# beta_t = s$beta_t
# test = s$test
# testn = t(apply(test, 1, function(x) x/sum(x)))
# 
# get_dist(wdfn[1,], beta_t)

mixedSCORE <- function(A, K, Thresh = NULL) {
  if(is.null(Thresh))
    Thresh = log(nrow(A))
  eig <- eigen(A) 
  lambdas <- eig$values[order(abs(eig$values), decreasing = T)[1:K]]
  Vs <- eig$vectors[,order(abs(eig$values), decreasing = T)[1:K]]
  R1 <- (Vs[,2:K,drop= F] / abs(Vs[,1]))
  Rhat <- sign(R1) * pmin(abs(R1), Thresh)
  
  # Vertex hunting
  require(combinat)
  Vlist = list()
  dList = list()
  i = 1
  for(L in (K+1):(3*K)) {
    km <- kmeans(Rhat, L)
    M <- km$centers
    combinations <- combn(1:L, K)
    dLs <- (apply(combinations, 2, function(x) {
      max(apply(M[-x,,drop = F],1, function(u) sqrt(get_dist(u, M[x,,drop = F]))))
    }))
    jmin <- which.min(dLs)
    Vlist[[i]] = M[combinations[,jmin], ,drop = F]
    dList[[i]] = min(dLs)
    i = i+1
  }
  deltaList = list()
  i  = 1
  for(L in (K+2):(3*K)) {
    V1 = Vlist[[i]]
    V2 = Vlist[[i+1]]
    deltaList[[i]] = min(sapply(permn(1:K), function(u) {
      max(sapply(1:K, function(v) sqrt(sum((V1[v,,drop = F]- (V2[u,,drop = F])[v,,drop = F])^2))))
    }))
    i = i+1
  }
  Ln = unlist(deltaList)/(1 + unlist(dList[2:length(dList)]))
  vk =  Vlist[[1 + which.min(Ln)]]
  # Membership reconstruction
  b1 = 1/sqrt(diag(vk %*% diag(lambdas[2:K,drop=F], nrow = K-1 ) %*% t(vk)) + lambdas[1])
  What = solve(rbind(t(vk), 1), rbind(t(Rhat), 1))
  pihatstar <- t(pmax(diag(1/b1) %*% What,0))
  pihat <- t(apply(pihatstar, 1, function(x) x/sum(x)))
  return(pihat)
}


# X <- Rhat
# chull(X)
# ## Not run: 
# # Example usage from graphics package
# plot(X, cex = 0.5)
# points(vk, col = "red")
# hpts <- chull(X)
# hpts <- c(hpts, hpts[1])
# lines(X[hpts, ])
# 

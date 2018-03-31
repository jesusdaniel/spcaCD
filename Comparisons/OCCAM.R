

occam <- function(A, K, Thresh = NULL) {
  E = eigen(A)
  Lhat = diag(E$values[1:K]*(E$values[1:K]>0))
  U = E$vectors[,1:K]
  ULU = U %*% Lhat %*% t(U) 
  ULhalf = U %*% sqrt(Lhat)
  alphahat = sum(A)/(nrow(A) * (nrow(A) - 1)*K)
  tau = 0.1 * alphahat^(0.2) * K^1.5/(nrow(A)^(0.3))
  X = t(apply(ULhalf,1, function(x) x/(sqrt(sum(x^2)+tau))))
  #install.packages("flexclust")
  require(flexclust)
  S = kcca(X, K,family = kccaFamily("kmedians"))
  #S = S@centers
  
  #Zhat = X%*%t(S)%*%solve((S%*%t(S)))
  #Zhat = t(apply(Zhat,1, function(x) x/(sqrt(sum(x^2)))))

  Zoccam = (X) %*% solve(parameters(S))
  Zoccam = Zoccam*(Zoccam>0)
  Zoccam = t(apply(Zoccam,1, function(x) x/(sum(x))))
  return(Zoccam)
}

distance_subs <- function(Z1, Z2) {
  base::norm(Z1 %*% solve(crossprod(Z1)) %*% t(Z1) - Z2 %*% solve(crossprod(Z2)) %*% t(Z2), "2")
}
add_noise <- function(Z, err=0.1) {
  Znoise <- err*matrix(runif(nrow(Z)*ncol(Z)), nrow = nrow(Z))
  t(apply(Z+Znoise, 1, function(x) x/sum(x)))
}

K= 3
N = 30
Z = matrix(0, N, K)

rho = 0.1
B = matrix(rho, K, K)
diag(B) = 1

Z[1:10,1] = 1
Z[11:20,2] = 1
Z[21:30,3] = 1

W = 0.5 * Z %*% B %*% t(Z)

Z0 = 0*Z
u1 = which(Z[,1] > 0)
u2 = which(Z[,2] > 0)
u3 = which(Z[,3] > 0)
Z0[u1, 1] = abs(eigen(W[u1, u1])$vectors[,1])
Z0[u2, 2] = abs(eigen(W[u2, u2])$vectors[,1])
Z0[u3, 3] = abs(eigen(W[u3, u3])$vectors[,1])

Z0 = Z
Zinit = add_noise(Z0,3)
eigen(crossprod(Zinit))
u1 =spca_iterative_DC(A = W, K = K, lambda = 0.7, lnorm = 2, 
                      thresh = hard_thresh,Z0 = Zinit, num_iter = 30,
                      TOL = 1e-10)
u1$iterations

u2 =spca_iterative_DC(A = W, K = K, 
                      lambda = 0.1, lnorm = 2, thresh = hard_thresh,
                      Z0 = add_noise(Z0,0.01), num_iter = 50, TOL = 1e-10)

u2$iterations
u2$Z

u1 = eigen(A[c(1:5, 11:20), c(1:5, 11:20)])
u1$vectors[,1]
u2$Z

Zhat = u2$Z
Z0
Zhat[,1]/Zhat[1,1]
Z

# now with degrees
theta = runif(N)
theta = rep(0.5,N)

Z = matrix(0, N, K)
Z[1:15,1] = 1
Z[11:20,2] = 1
Z[21:30,3] = 1
Z = t(apply(Z, 1, function(x) x/sum(x)))

W = diag(theta) %*% Z %*% B %*% t(Z) %*% diag(theta)

Z0 = 0*Z
u1 = which(Z[,1] > 0)
u2 = which(Z[,2] > 0)
u3 = which(Z[,3] > 0)
Z0[u1, 1] = abs(eigen(W[u1, u1])$vectors[,1])
Z0[u2, 2] = abs(eigen(W[u2, u2])$vectors[,1])
Z0[u3, 3] = abs(eigen(W[u3, u3])$vectors[,1])
B0 = solve(crossprod(Z0)) %*% crossprod(Z0, W)%*% Z0 %*% solve(crossprod(Z0)) 
Z0 %*% diag(sqrt(diag(B0)))
Z0 %*% solve(B0)
distance_subs(Z0, Z)

distance_subs(Z0, eigen(W)$vectors[,1:K])

Znoise<- add_noise(Z0, err = 0.3)
Zspca = spca_iterative_DC(A = W, K = K, lambda = 0.15, 
                          lnorm = 2, thresh = soft_thresh,Z0 = Z, num_iter = 30)

Zspca$Z
distance_subs(Zspca$Z, eigen(W)$vectors[,1:K])

Znoise <- add_noise(Z0, err = 5)
Z0 <- add_noise(Z0, err = 0.3)
spca_iterative(A = W, K = K, lambda = 0.6, lnorm = 1, Z0 = Z, num_iter = 30)
spca_iterative(A = W, K = K, lambda = 0.7, lnorm = 1, Z0 = Z, num_iter = 30, thresh = hard_thresh)
Zinverse <- spca_iterative_inverse(A = W, K = K, lambda = 0.5, lnorm = 1, Z0 = Z0, num_iter = 30, 
                       TOL = 1e-09,
                       thresh = hard_thresh)$Z
Zinverse
distance_subs(Z0, eigen(W)$vectors[,1:K])
distance_subs(Zinverse, eigen(W)$vectors[,1:K])

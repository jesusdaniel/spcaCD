## 
p = 0.3
q = 0.15
Psi = (p-q)*diag(K)  + q
nk = 150
N = K*nk
A = kronecker(Psi,matrix(rep(1,(nk)^2),nrow = nk))
A[101:200, 1:100] = 0.8 * 0.3 + 0.2 * 0.15
A[1:100, 101:200] = 0.8 * 0.3 + 0.2 * 0.15
A[101:200, 101:200] = 0.8 * (0.8 * 0.3 + 0.2 * 0.15) +
  0.2 * (0.2 * 0.3 + 0.8 * 0.15)
A[101:200, 201:300] = 0.2 * 0.3 + 0.8 * 0.15
A[201:300, 101:200] = 0.2 * 0.3 + 0.8 * 0.15

A[101:200,] = 0.18
A[,101:200] = 0.18
levelplot(A)

# SBM

K = 3
N = 500
A <- sbm(0.5, 0.1, K, N/K)
Ztrue = A$membership
#levelplot(A$Aobs)
A = A$Aobs
A <- generate_overlapping(n = 500)
Ztrue <- A$membership
A <- A$Aobs
Z0 <- initialization(N, K)

Z


eig <- eigen(A)$vectors[,1:K]
lambdas <- eigen(A)$values
v1 <- eigen(A[1:(N/K), 1:(N/K)])$vectors[,1]
v2 <- abs(eigen(A[(N/K+1):N,(N/K+1):N])$vectors[,1])

Vs = cbind(c(v1,rep(0, N/K)), c(rep(0,N/K), v2))
crossprod(Vs)


table(binarize(Z)[,1],Ztrue)


# rotated eigenvectors
theta = 3*pi/4
Rot = matrix(c(cos(theta), sin(theta), - sin(theta), cos(theta)), 2)
plot(eig[,1], eig[,2], col = Ztrue, xlim = c(-0.1, 0.1), ylim = c(-0.1, 0.1))
abline(h=0); abline(v=0)
Rot = matrix(c(cos(theta), sin(theta), - sin(theta), cos(theta)), 2)
eigr =  (eig %*% diag(c(3,1))) %*% Rot
plot(eigr[,1], eigr[,2], col = Ztrue, xlim = c(-0.1, 0.1), ylim = c(-0.1, 0.1))
abline(h=0); abline(v=0)


Z0 =  (eig %*% diag(c(3,1))) %*% Rot






# rank 0
sbm_B <- function(B,  nk) {
  Psi = B
  A = kronecker(Psi,matrix(rep(1,(nk)^2),nrow = nk))
  Aobs = apply(A,MARGIN = c(1,2),function(u) rbinom(1,1,prob = u))  
  Aobs[lower.tri(Aobs)] = 0
  Aobs =  Aobs + t(Aobs)
  diag(Aobs) = 0
  return(list(Aobs= Aobs, membership = kronecker(1:K,rep(1,nk)), W = A))
}

B = matrix(0, 3,3)
diag(B) = 1
ps = 1/sqrt(2) #+ 0.1
B[2:3,1] = B[1, 2:3] = ps
B[2,3] = B[3,2] = 0#2*ps^2-1
B
2.56/4
1.28/2

0.8/0.64
eigen(B)
B= 0.5*B
A = sbm_B(B, 200)
Vs = eigen(A$Aobs)$vectors
(eigen(A$Aobs)$values)[1:5]
plot(Vs[,1])
plot(Vs[,2])
plot(Vs[,3])
plot(Vs[,600])

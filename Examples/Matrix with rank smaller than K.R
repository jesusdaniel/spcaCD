D = matrix(c(1,9,9,1), 2)

D = matrix(0.1, 3, 3)
diag(D) = 1
D
w1 = c(0.5,0.5,0)
a = solve(D) %*% w1
s1 <- crossprod(a,w1) 

w1
D %*% a

B = rbind( cbind(D, w1/sqrt(s1)), c(w1 / sqrt(s1),1))
#B = rbind( cbind(D, w1), c(w1, s1))

rankMatrix(B)


D = B[2:4, 2:4]
w1 = D[1:3,1]
a = solve(D) %*% w1
s1 <- crossprod(a,w1) 


B
u = lm(B[,4] ~ B[,1:3] - 1)
u = lm(B[,1] ~ B[,2:4] - 1)
u = lm(B[,2] ~ B[,c(1,3:4)] - 1)
u = lm(B[,3] ~ B[,c(1:2,4)] - 1)
u$residuals
u$coefficients
eigen(B)
B

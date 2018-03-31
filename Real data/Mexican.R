## SBM

library(igraph)

#DipNet = read.graph(file = "C:\\Users\\JDAR\\Box Sync\\Etc\\Politica\\Diputados_und.graphml",format = "graphml")

A = A_diputados
A <- A_diputados * t(A_diputados)
A[1:5,1:5]
isSymmetric(A)
degrees1 = apply(A,1,sum)
degrees2 = apply(A,2,sum)

nonzero = which(degrees1 >0 | degrees2 >0)
A = A[nonzero,nonzero]
A2 = A%*%A
nonzero2 = which(apply(A2,1,sum)!=1)
A = A[nonzero2,nonzero2]
dim(A)

diag(A) = 0
sum(diag(as.matrix(A)))
A[which(as.matrix(A)>1,arr.ind=T)]=1
sum(A)/2
dim(A)

## Communities
party <- ifelse(is.na(as.character(diputados$Partido)), "NA", 
                as.character(diputados$Partido))
names <- df_diputados$name
names <- (names[nonzero])[nonzero2]
party_Fac = as.character((party[nonzero])[nonzero2])
party <- as.numeric(factor(party_Fac))

V = eigen(A)
plot(V$values)
plot(V$vectors[,1],V$vectors[,2],col = party+1)
plot(V$vectors[,2]/V$vectors[,1],col = party+1)

degree = apply(A,1,sum)
mean(degree)
plot(degree)
hist(degree,breaks = 50)

source("community_detection.R")
source("sparsePCA.R")
sPCA = overlapping_community_detection(A = A,K=10,0.9,iterations = 500)
plot(sPCA$Zraw,col = party+1)

table(sPCA$membership %*% (1:10), party_Fac)
sPCA$membership

Lap = laplacian(A,tau = 0.00)
sPCA = overlapping_community_detection(A = A, K=10,0.9,iterations = 500)
community = apply(sPCA$Zest,1,which.max)
table(community,party_Fac)
#with tau = 0, ...
#with tau = 0.001, best error is 59 (lambda around 0.00441)
#with tau=0.005, best error is 68 (with 0.0074 lambda)
#with tau = 0.007, best error is 63 (lambda around 0.0095)
#with tau = 0.01, best error is 65 (lambda around 0.0105)

plot(sPCA$Zraw,col = party+1)

overlapping_nodes <- which(apply(sPCA$membership,1, sum) >1)
names[overlapping_nodes]

lambdapath = seq(0.2,0.95, 0.05)

result = sparsePCA_comdet_path(A, K, lambdapath, DC = F, Ztrue = Z,
                               l_norm = 1, num_trials = 5, maxiter = 100,TOL=1e-4)

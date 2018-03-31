K = 2
rho = 0.1
  B = array(rho,dim = c(K,K))
  diag(B) = 1 
npure = 5
noverlap = 10
Z = Matrix(0,nrow = K*npure + noverlap, ncol=K)
Z[1:npure,1] = 1
Z[(npure+1):(2*npure),2] = 1

x = seq(0.1,0.9, length.out = noverlap)
Z[(K*npure + 1):(K*npure + noverlap),c(1,2)] = cbind(x, 1-x)
d= 5

A = Z%*%B%*%t(Z)
n = nrow(A)
alpha = sum(apply(A,2,sum))/(n*d)
W = A/alpha
W = pmin(W,1)
Aobs = apply(W,MARGIN = c(1,2),function(u) rbinom(1,1,prob = u))
Aobs[lower.tri(Aobs)] = 0
Aobs =  Aobs + t(Aobs)
diag(Aobs) = 0
Ztrue = Z
A = Aobs
A = W

# score---------------------------------------
sol_mSCORE <- mixedSCORE(Aobs, 2)
plot(sol_mSCORE[1:(K*npure),])
points(sol_mSCORE[(K*npure+1):nrow(A),], col = "red")
bestSCORE = best_exNIV_threshold(sol_mSCORE, Z)
exNIV_threshold(sol_mSCORE, Z)
accuracy(binarize(sol_mSCORE), binarize(Z))
accuracy(binarize(sol_mSCORE[1:(K*npure),]), binarize(Z[1:(K*npure),]))
bestSCORE
#occam----------------------------------------
sol_occam <- occam(A, K)
plot(sol_occam[1:(K*npure),])
points(sol_occam[(K*npure+1):nrow(A),], col = "red")
bestOCCAM = best_exNIV_threshold(sol_occam, Z)
exNIV_threshold(sol_occam, Z)
bestOCCAM
accuracy(binarize(sol_occam), binarize(Z))
accuracy(binarize(sol_occam[1:(K*npure),]), binarize(Z[1:(K*npure),]))

# SPCA----------------------------------------
lambdapath = seq(0.3, 0.99, length.out = 20)
sol_spca <- sparsePCA_comdet_path(A = A, K = K, lambdapath = lambdapath, DC = F, Ztrue = Z, 
                      l_norm = 1, num_trials = 5)
sol_spca$accuracy_estimated
sol_spca$exNVIs_estimated
sol_spca$best_exNVI
sol_spca$best_accuracy
solspca = sol_spca$Z_hat_SPCA
plot(solspca[1:(K*npure),])
points(solspca[(K*npure+1):nrow(A),], col = "red")
accuracy(binarize(solspca[1:(K*npure),]), binarize(Z[1:(K*npure),]))
accuracy(binarize(solspca), binarize(Z))

bestSPCA = 0
for(lambda in seq(0.35, 0.99, 0.01)) {
  Z0 <- initialization(n, K)
  lambda = 0.7
  sol_spca <- spca_iterative(A, K, lambda = lambda, lnorm = 1, Z0 = Z0, 
                                num_iter = 50, TOL = 1e-08, save_path = T)
  bestSPCA = max(bestSPCA, exNiV((sol_spca$Z>0), Z>0))
  solspca <- sol_spca$Z
  print(accuracy(binarize(solspca), binarize(Z)))
  #print(exNiV((sol_spca$Z>0), Z>0))
}


sol_spca <- spca_iterative_DC3(A, K, lambda = 0.05, lnorm = 2,
                               Z0 = add_noise(Ztrue), thresh = hard_thresh,
                               num_iter = 50, TOL = 1e-06, save_path = T)

sol_spca <- spca_iterative_DC3(A, K, lambda = 0.2, lnorm = 2,
                               Z0 = initialization(20,2), thresh = hard_thresh,
                               num_iter = 50, TOL = 1e-06, save_path = T)

sol_spca$Z
sol_spca$iterations
distance_subs(sol_spca$Z, Ztrue)
eigen(A[c(1:5, 18:20), c(1:5, 18:20)])$vectors[,1]

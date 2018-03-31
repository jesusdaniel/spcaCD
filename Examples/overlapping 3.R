K = 3
N = 500
npure = 150
noverlap = 1
N = K*(npure + noverlap)
N = 350 + 300 + K*noverlap
d = 40
rho = 0.3
B = array(rho,dim = c(K,K))
diag(B) = 1 
B[1,2] = B[2,1] = 0.3
A <- generate_overlapping_overlapsize(rho = rho, d = d, B = B,
                                      npure = npure, noverlap = noverlap)
A <- generate_overlapping_comsize(rho = rho, d = d, 
                                      npure1 = 350, noverlap = noverlap)
#A <- generate_overlapping(n = 500, rho = rho, d = d)
Ztrue <- A$membership
A <- A$Aobs
N = nrow(A)
Z0 <- initialization(N, K)
colSums(Ztrue)

# bestSPCAinv = 0
# for(lambda in seq(0.8, 0.99, 0.01)) {
#   Z0 <- initialization(N, K)
#   sol_spca_inv <- spca_iterative_inverse(A, K, lambda = 0.5, 
#                                          lnorm = 1, 
#                              Z0=Z0, num_iter = 50, TOL = 1e-06, save_path = T)
#   print(exNiV((sol_spca_inv$Z>0), Ztrue>0))
#   bestSPCAinv = max(bestSPCAinv, exNiV((sol_spca_inv$Z>0), Ztrue>0))
# }

bestSPCA = 0
for(lambda in seq(0.1, 0.99, 0.01)) {
  Z0 <- initialization(N, K)
  sol_spca <- spca_iterative(A, K, lambda = lambda, lnorm = 1, 
                              Z0=Z0, num_iter = 50, TOL = 1e-04, save_path = T)
  print(exNiV((sol_spca$Z>0), Ztrue>0))
  bestSPCA = max(bestSPCA, exNiV((sol_spca$Z>0), Ztrue>0))
}

bestSPCA_DC = 0
lambda=0
Z0= Ztrue
spcadc_exNIV_path = c()
lambda_path <- seq(0.01, 0.99, 0.01)
Zpath = list()
i=1
for(lambda in lambda_path) {
  Z0 <- initialization(N, K)
  #Z0= sol_mSCORE
  sol_spca <- spca_iterative_DC(A, K, lambda = lambda, lnorm = 2,
                                Z0 = Z0, thresh = hard_thresh,
                                num_iter = 50, TOL = 1e-04, save_path = T)
  
  print(exNiV((sol_spca$Z!=0), Ztrue>0))
  #print(distance_subs(eigen(A)$vectors[,1:K], sol_spca$Z))
  Zpath[[i]] = sol_spca$Z
  i=i+1
  spcadc_exNIV_path = c(spcadc_exNIV_path,exNiV((sol_spca$Z!=0), Ztrue>0))
  
}
plot(spcadc_exNIV_path)
max(spcadc_exNIV_path)
sol_spca$Z
distance_subs(eigen(A)$vectors[,1:K], sol_spca$Z)
distance_subs(eigen(A)$vectors[,1:K], Ztrue)

bestSPCA_DC_soft = 0
for(lambda in seq(0.35, 0.99, 0.01)) {
  Z0 <- initialization(N, K)
  sol_spca <- spca_iterative_DC(A, K, lambda = lambda, lnorm = 2, thresh = soft_thresh,
                                Z0=Z0, num_iter = 50, TOL = 1e-04, save_path = T)
  print(exNiV((sol_spca$Z>0), Ztrue>0))
  bestSPCA_DC_soft = max(bestSPCA_DC_soft, exNiV((sol_spca$Z>0), Ztrue>0))
}

############################################
sol_mSCORE <- mixedSCORE(A, K)
bestSCORE = best_exNIV_threshold(sol_mSCORE, Ztrue)
exNIV_threshold(sol_mSCORE, Ztrue)
accuracy(binarize(sol_mSCORE), binarize(Ztrue))
accuracy(binarize(sol_mSCORE[1:(K*npure),]), binarize(Ztrue[1:(K*npure),]))
bestSCORE
#occam----------------------------------------
sol_occam <- occam(A, K)
plot(sol_occam[,1])
plot(sol_occam[,2])
plot(sol_occam[,3])
bestOCCAM = best_exNIV_threshold(sol_occam, Ztrue)
exNIV_threshold(sol_occam, Ztrue)
bestOCCAM
accuracy(binarize(sol_occam), binarize(Ztrue))
accuracy(binarize(sol_occam[1:(K*npure),]), binarize(Z[1:(K*npure),]))

bestSCORE
bestOCCAM
bestSPCA
bestSPCA_DC
bestSPCA_DC_soft

source("./SPCA_Rcode/community_detection.R")
source("./SPCA_Rcode/evaluate_estimation.R")
source("./SPCA_Rcode/likelihood_functions.R")
source("./SPCA_Rcode/sparsePCA.R")
source("./SPCA_Rcode/sparsePCA_comdet_path.R")

args = commandArgs(TRUE)

A = as.matrix(read.csv(sprintf('./SPCA_temp/temp_A.csv'),
                       header=FALSE))


Z = as.matrix(read.csv(sprintf('./SPCA_temp/temp_Ztrue.csv'),
                        header=FALSE))

K = ncol(Z)

lambdapath = seq(0.2,0.95, 0.05)

result = sparsePCA_comdet_path(A, K, lambdapath, DC = F, Ztrue = Z,
                               l_norm = 1, num_trials = 5, maxiter = 100,TOL=1e-4)

result_DC = sparsePCA_comdet_path(A, K, lambdapath, DC = T, Ztrue = Z,
                               l_norm = 1, num_trials = 5, maxiter = 100,TOL=1e-4)

write.table(as.matrix(result$Z_hat_SPCA), file=sprintf('./SPCA_temp/Z_hat_SPCA.csv'),
            sep=',', row.names=FALSE, col.names=FALSE)

write.table(as.matrix(result_DC$Z_hat_SPCA_DC), file=sprintf('./SPCA_temp/Z_hat_SPCA_DC_BICDC.csv'),
            sep=',', row.names=FALSE, col.names=FALSE)
write.table(as.matrix(result_DC$Z_hat_SPCA), file=sprintf('./SPCA_temp/Z_hat_SPCA_DC.csv'),
            sep=',', row.names=FALSE, col.names=FALSE)

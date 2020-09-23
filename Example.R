# Load R functions and packages
source("R/loadAll.R")


# Load Mexican representatives data ---------------------
load("data/MexicanRepresentatives.RData")
library(lattice)

# Node attributes
name_representatives
colors_partidos
party_representatives
unique(party_representatives)

# Plot adjacency
levelplot(Adjacency_representatives)

# Generate common layout to reuse later in all plots
net <- graph.adjacency(Adjacency_representatives, mode = "undirected", diag = F)
layout.fr <- layout.fruchterman.reingold(net)

# Wrappper of plot.igraph function to color vertices by community membership
membership_plots(A = Adjacency_representatives, 
                communities = factor(party_representatives), 
                colors = colors_partidos,
                nodesizes = rep(4, n), 
                layoutpos = layout.fr)





# Fit overlapping communities ------------------------------------------

# Make undirected network by symmetrizing adjacency
A  <- Adjacency_representatives
A <- 1*(A + t(A) > 0 )

n <- ncol(A) #number of nodes



# Run Algorithm 1 (method = "SPCAeig") overlapping communities with degree heterogeneity
comdet_alg1 <- overlapping_community_detection(A = A, 
                                               K = 9, # chosen as the number of parties
                                               method = "SPCAeig", 
                                               lambda = 0.9, 
                                               initial = 1,   # Random initialization
                                               num_trials = 10) # number of different random initializations

membership_pie_plots(A = A, 
                     Z = comdet_alg1$Zest,  # membership matrix
                     nodesizes = rep(4, n), 
                     layoutpos = layout.fr)


# Run Algorithm 2 (method = "SPCA_Z") overlapping communities with homgeneous degrees
comdet_alg2 <- overlapping_community_detection(A = A, 
                                               K = 9, # chosen as the number of parties
                                               method = "SPCA_Z", 
                                               lambda = 0.9, 
                                               initial = 1,   # Random initialization
                                               num_trials = 10) # number of different random initializations

membership_pie_plots(A = A, 
                     Z = comdet_alg2$Zest,  # membership matrix
                     nodesizes = rep(4, n), 
                     layoutpos = layout.fr)



# Run Algorithm 2 with 3 communities for illustration
#  initial = 3 uses SCORE (non-overlapping community detection) for initialization
comdet_alg2_K3 <- overlapping_community_detection(A = A, 
                                               K = 3, # chosen as the number of parties
                                               method = "SPCA_Z", 
                                               lambda = 0.9, 
                                               initial = 3) # number of different random initializations
membership_pie_plots(A = A, 
                     Z = comdet_alg2_K3$Zest,  # membership matrix
                     nodesizes = rep(4, n), 
                     layoutpos = layout.fr)


# Plot membership paths
# Note: In order to get consistent memberships for different lambda values, "initial" should be
# equal to 2 (user-specified initial value) and provide initial matrix "Z0", or 
# equal to 3 (non-overlapping communities obtained by SCORE)

# Sequence of thresholding values
lambdapath <- seq(0.99, 0.2, by = -0.01)

comdet_path2_K3 <- SPCA_CD_path(A = A, 
                                K = 3,
                                method = "SPCA_Z", 
                                lambdapath = lambdapath,
                                thresh = soft_thresh, # the soft thresholding function creates smoother node paths
                                initial = 3)

# The following plot shows the node membership paths (see Figure 5.3 in https://arxiv.org/pdf/2009.10641.pdf)
# Each line corresponds to a different node. 
# A non-zero value indicates node membership to the community in the corresponding plot 
plot_nodepaths(comdet_path = comdet_path2_K3, 
               lambdas = lambdapath, 
               colors = colors_partidos[factor(party_representatives)])

# Select best model by BIC
cd.bic <- SPCA_CD_path.BICselection(A,comdet_path2_K3)
Z_alg2_bic <- cd.bic$Z_minBIC$Zest
cd.bic$lambda_minBIC
membership_pie_plots(A = A, 
                     Z = Z_alg2_bic,
                     nodesizes = rep(4, n), 
                     layoutpos = layout.fr)


# Select best model by network cross-validation
cv.path_eig <- cv.SPCA_CD_path(A = A, 
                               K = 3, 
                               method = "SPCA_Z", 
                               impute_function = 3,
                               loss_function = 2,
                               lambdapath = lambdapath,
                               initial = 3)

Z.cv = cv.path_eig$Z_minCV$Zest
membership_pie_plots(A = A, 
                     Z = Z.cv,
                     nodesizes = rep(4, n), 
                     layoutpos = layout.fr)

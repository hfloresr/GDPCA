#!/usr/bin/env Rscript

.libPaths("/project/k01/markomg/R/x86_64-pc-linux-gnu-library/3.3/")

library(iterators)
library(doParallel)
library(R.matlab)
library(gdpc)
library(MASS)


# Read in pre_pmcao data
data = readMat('data/F141020-lfp-5min-1kHz.mat')
Z = data$pre.pmcao


# Channels to exclude from Z
excludes = c(1:5, 7:8, 11:12, 14:16, 18:26, 28:32)

ncores = 30
cl <- makeCluster(ncores)
registerDoParallel(cl)

clusterEvalQ(cl, .libPaths("/project/k01/markomg/R/x86_64-pc-linux-gnu-library/3.3/"))

k=20  # lag
num_epochs = 300  # All time points
epoch_models = foreach (i=1:num_epochs, .packages='gdpc', .combine='cbind') %dopar%{
        model <- gdpc(Z[((i-1)*1000+1):(i*1000), excludes], k)
        model
}
stopCluster(cl)
dump("epoch_models", sprintf("pre_epoch_models%d.Rdmpd", k))

# Write matrix
#colnames(epoch_factors) = NULL
#write.matrix(epoch_factors,'gdpca_factors_k10.rmat', sep='\t')

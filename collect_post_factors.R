#!/usr/bin/env Rscript

.libPaths("/project/k1258/floresh/R/x86_64-pc-linux-gnu-library/3.3/")

library(doParallel)
library(R.matlab)
library(gdpc)
#library(MASS)

# Read in pre_pmcao data
data = readMat('data/F141020-lfp-5min-1kHz.mat')
Z = data$post.pmcao

# Channels to exclude from Z
excludes = c(1:5, 7:8, 11:12, 14:16, 18:26, 28:32)
exc_chs_11_15_16 = c(54, 63, 64, 114:117, 136:141, 151:154, 161:168,
								 182:187, 200:208, 223:225, 237:240, 296:297)
exc_chs_15_16 = c(62, 65, 66)
exc_chs_11_16 = c(236, 241, 243)
exc_chs_16 = c(171, 172)

ncores = 32
cl <- makeCluster(ncores)
registerDoParallel(cl)

clusterEvalQ(cl, .libPaths("/project/k1258/floresh/R/x86_64-pc-linux-gnu-library/3.3/"))

k=20  # lag
num_epochs = 300  # All time points
epoch_models = foreach (i=1:num_epochs, .packages='gdpc', .combine='cbind') %dopar%{
    if(i %in% exc_chs_11_15_16) {
		    exs = exc_chs_11_15_16
    } else if( i %in% exc_chs_15_16) {
				exs = exc_chs_15_16
    } else if(i %in% exc_chs_11_16) {
				exs = exc_chs_11_16
    } else if(i %in% exc_chs_16) {
				exs = exc_chs_16
    } else{
			exs = excludes
    }
    model <- gdpc(Z[((i-1)*1000+1):(i*1000), exs], k)
    model
}
stopCluster(cl)
dump("epoch_models", sprintf("post_epoch_models.Rdmpd", k))

library(doParallel)
library(R.matlab)
library(GDPCA)
library(MASS)

# Read in pre_pmcao data
data = readMat('data/F141020-lfp-5min-1kHz.mat')
pre_pmcao = data$pre.pmcao

# Get factors for all 300 epochs for pre_pmcao
num_epochs = 300
k = 10  # lags
excludes = c(1:5, 7:8, 11:12, 14:16, 18:26, 28:32)

#ncores = 34
ncores = 2
cl <- makeCluster(ncores)
registerDoParallel(cl)

epoch_factors = foreach (i=1:num_epochs, .packages='gdpc', .combine='cbind') %dopar%{
  fit <- gdpc(pre_pmcao[((i-1)*1000+1):(i*1000), excludes], k)
  factor = (fit$f - mean(fit$f)) / sd(fit$f)
  factor
}
stopCluster(cl)

# Write matrix
colnames(test_epoch_factors) = NULL
write.matrix(epoch_factors,'gdpca_factors_k10.rmat', sep='\t')

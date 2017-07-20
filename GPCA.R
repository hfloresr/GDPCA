# This is a simple implementation of the dynamic principal component anaylsis
# methods presented in the following paper.
#
# Pe?a, Daniel, and Victor J. Yohai. "Generalized Dynamic Principal Components."
# Journal of the American Statistical Association 111.515 (2016): 1121-1131.
#
# TODOs:
# - (1.) vectorize some operations since current implementation using
# may explicit for loops which can be extremely inefficient
# - (2.) refactor the code and wrap them in a class for this model.
# The class may contains methods like model.train(), model.evaluate(), ...
# - (3.) add test cases (currently zero test case)

library(R.matlab)

init_f <- function(Z, k){
  n = ncol(Z)
  f_init = matrix(0, n + k, 1)
  Z_t = t(Z)

  mzt = rowMeans(Z_t)
  mzt = as.matrix(mzt)

  Z_new = matrix(0, dim(Z_t)[1], dim(Z_t)[2])

  for(i in 1:dim(z)[2]) {
    Z_new[, i] = Z_t[, i] - mzt
  }

  #z_new = sweep(z, 2, rowMeans(z), "-")
  Z_svd = svd(z_new)
  f_init[1:n, ] = Z_new %*% as.matrix(Z_svd$v[, 1])

  f_init = f_init/sd(f_init)
  return (f_init)
}


run_gdpc <- function(Z, k, f_ini=NULL, tol=1e-04, maxiter=100) {
  # Each ROW is a DIFFERENT time series
  # run the DPC model
  m = nrow(z);
  T = ncol(z);


  f = ifelse(is.null(f_ini), init_f(Z, k), f_ini)
  alpha = NaN
  beta = NaN
  loss_values = matrix(0, nrow=maxiter, ncol=1)

  for (train_iter in 1:maxiter) {
    a = run_train_step(z, k, f, alpha, beta);
    f = a$f_new;
    alpha = a$alpha_new;
    beta = a$beta_new
    loss_values[train_iter] = evaluate_op(z, k, f, alpha, beta);
  }

  ## plot result
  #plot(loss_values, lwd = 1.5, main="train loss (reconstruction MSE)", xlab = "iteration");
  #print(loss_values)

  #plot(f, lwd = 1.5, main = "Factor", xlab = "time" );
  #lines(f, lty = 1)
}

# ---------- below are all helper functions -------------------------------
# not optimized for speed yet (use a lot of for loops)
#

run_train_step <- function(z, k, f, alpha, beta) {
  m = nrow(z);
  T = ncol(z);
  ab = alpha_beta(z, f, k);
  alpha_new = ab$alpha;
  beta_new = ab$beta;
  f_star = f_alpha_beta(z, k, alpha_new, beta_new);
  f_centered = f_star -  mean(f_star);
  f_new = sqrt(T + k) * f_centered / max(svd(f_centered)$d); ########
  return(list("f_new"=f_new, "alpha_new" = alpha_new, "beta_new" = beta_new))
}


evaluate_op <- function(z, k, f, alpha, beta){
  loss = mean_squared_error(z, k, f, alpha, beta)
  return (loss)
}


C_alpha <- function(z_j, alpha_j, k){
  T = ncol(z_j);
  C = matrix(0, T + k, k + 1);
  for (t in 1:(T + k)){
    for (q in 1:(k + 1)){
      if ((q >= max(1, t - T + 1)) && (q <= min(k + 1, t))){
        C[t, q] = z[t - q + 1] - alpha_j
      }
    }
  }
  return (C)
}


D_beta <- function(beta_j, T, k) {
  # Construct matrix D w.r.t beta_j, T, and k
  # beta_j is a row vector of the Beta matrix
  D = matrix(0, T+k, T+k)
  for (t in 1:T+k) {
    for (r in max(t-k, 1):min(t,T)) {
      for(q in r:r+k) {
        D[t,q] = D[t,q] + beta_j[q-r+1] %*% beta_j[t-r+1]
      }
    }
  }
  return (D)
}


f_alpha_beta <- function(Z, alpha, beta, k){
  # Get optimal f wrt Z, beta, alpha, and k
  m = nrow(Z)
  T = ncol(Z)
  D = matrix(0, T+k, T+k)
  f = matrix(0, T+k, 1)
  for (j in 1:m) {
    D = D + D_beta(beta[j,], T, k)
    f = f + C_alpha(Z[j,], alpha[j]) %*% beta[j,]
  }
  f <- ifelse(rcond(D)>1e-10, solve(D,f), pinv(D))
  return (f)
}


F_f <- function(Z, f, k){
  # Gets F matrix s.t. F(f) is T x (k+2) w/ t-th row (f_t, f_t+1,...,f_t+k, 1)
  m = nrow(Z);
  T = ncol(Z);
  F = matrix(0,nrow=T,ncol=k+2)
  for (t in 1:nrow(F)){
    F[t, ] = cbind(Conj(t(f[t:(t + k)])), 1);
  }
  return(F);
}


alpha_beta <- function(Z, f, k){
  # Find the optimal beta
  # input: Z data matrix, f-principal components, k-leads
  # output: beta matrix of loadings, alpha-intercepts
  m=nrow(Z)
  T=ncol(Z)
  F=F_f(Z, f, k)
  FtF_inv = solve(Conj(t(F)) %*% F);
  FtF_inv_Ft = FtF_inv %*% Conj(t(F));
  tmp = z %*% Conj(t(FtF_inv_Ft));
  alpha = Conj(t(tmp[, k + 2]));
  beta = tmp[, 1:(k + 1)];
  #res = Z - beta %*% F
  return(list("alpha"=alpha, "beta"=beta))
}


mean_squared_error <- function(z, k, f, alpha, beta){
  m = nrow(z);
  T = ncol(z);
  sum_squared_error = 0;
  for (j in 1:m){
    for (t in 1:T){
      z_jt_predict = alpha[j];
      for (i in 0:k){
        z_jt_predict = z_jt_predict + beta[j, i + 1] * f[t + i];
      }
      sum_squared_error = sum_squared_error + (z[j, t] - z_jt_predict)^2;
    }
  }
  err = sum_squared_error / (T * m);
  return(err);
}



#z = matrix(rnorm(10*100,mean=0,sd=1), 10, 100)
#zdata = readMat("z.mat")
Z = readMat('F141020-lfp-5min-1kHz.mat')
pre_pmcao = z$pre.pmcao
Z = pre_pmcao[((1-1)*1000+1):(1*1000), ]
Z_t = t(z)
k=2

run_gdpc(Z_t, k, f_ini=NULL, tol=1e-04, maxiter=1000)
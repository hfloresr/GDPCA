# This is a simple implementation of the dynamic principal component anaylsis
# methods presented in the following paper.
#
# Peña, Daniel, and Victor J. Yohai. "Generalized Dynamic Principal Components."
# Journal of the American Statistical Association 111.515 (2016): 1121-1131.
#
# TODOs: 
# - (1.) vectorize some operations since current implementation using
# may explicit for loops which can be extremely inefficient
# - (2.) refactor the code and wrap them in a class for this model. 
# The class may contains methods like model.train(), model.evaluate(), ...
# - (3.) add test cases (currently zero test case)

init_f <- function(z, k){
  
  n = ncol(z)
  f_init = matrix(0, n + k, 1)
  z = t(z)
  
  mzt = rowMeans(z)
  mzt = as.matrix(mzt)
  
  z_new = matrix(0, dim(z)[1], dim(z)[2])
  
  for(i in 1:dim(z)[2]){
    
    z_new[, i] = z[, i] - mzt
    
  }
  
  #z_new = sweep(z, 2, rowMeans(z), "-")
  z_svd = svd(z_new)
  f_init[1:n, ] = z_new %*% as.matrix(z_svd$v[, 1])
  #d_m = vec2diag(z_svd$d)
  #temp_m = z_svd$u %*% d_m
  #temp_m = as.matrix(temp_m[, 1])
  
  #for(i in 1:n) {
  
  # f_init[i,1] = temp_m[i, 1]
  
  #}
  
  f_init = f_init/sd(f_init)
  
  return(f_init)
}

dynamic_principal_components <- function() {
  # run the DPC model
  
  # set up
  
  #z = matrix(rnorm(10*100,mean=0,sd=1), 10, 100) 
  #zdata = readMat("z.mat")
  z = pre_pmcao[((1-1)*1000+1):(1*1000), ]
  z = t(z)
  k1 = 50;
  k2 = 30;
  k = k1 + k2;
  m = nrow(z);
  T = ncol(z);
  #f_init = matrix(rnorm((T+k)*1,mean=0,sd=1), T + k, 1)
  #f_initdata = readMat("f.mat")
  f_init = matrix(0, T + k, 1)
  alpha_init = NaN;
  beta_init = NaN;
  
  # train
  #f = init_f(z, k);
  f = f_init;
  alpha = alpha_init;
  beta = beta_init;
  train_iterations = 20;
  loss_values = matrix(0, nrow=train_iterations, ncol=1)
  
  
  for (train_iter in 1:train_iterations) {
    a = run_train_step(z, k, f, alpha, beta);
    f = a$f_new;
    alpha = a$alpha_new;
    beta = a$beta_new
    loss_values[train_iter] = evaluate_op(z, k, f, alpha, beta);
  }
  
  # plot result
  
  plot(loss_values, lwd = 1.5, main="train loss (reconstruction MSE)", xlab = "iteration");
  print(loss_values)
  
  
  plot(f, lwd = 1.5, main = "Factor", xlab = "time" );
  lines(f, lty = 1)
}

# ---------- below are all helper functions -------------------------------
# not optimized for speed yet (use a lot of for loops)
#

run_train_step <- function(z, k, f, alpha, beta){
  
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


evaluate_op <-function(z, k, f, alpha, beta){
  loss = mean_squared_error(z, k, f, alpha, beta)
  return(loss)
}


C_alpha <- function(z, alpha, k){
  m = nrow(z);
  T = ncol(z);
  C = array(0, dim=c(m, T + k, k + 1));
  for (j in 1:m){
    for (t in 1:dim(C)[2]){
      for (q in 1:dim(C)[3]){
        if ((q >= max(1, t - T + 1)) && (q <= min(k + 1, t))){
          C[j, t, q] = z[j, t - q + 1] - alpha[j]; #####
        }
      }
    }
  }
  return(C)
}

D_beta = function(z, beta, k){
  m = dim(z)[1]
  T = dim(z)[2]
  D = array(0,dim=c(m,T + k,T + k))
  for (j in 1:m) {
    for (t in 1:dim(D)[2]){
      for (q in 1:dim(D)[3]){
        if ((q >= max(t - k, 1)) && (q <= min(t + k, T + k))) {
          for (v in max(1, t - k, q - k):min(t, q, T)){
            D[j, t, q] = D[j, t, q] + beta[j, q - v + 1] * beta[j, t - v + 1]; 
          }
        }
      }
    }
  }
  
  return(D)
}


f_alpha_beta <- function(z, k, alpha, beta){
  m = nrow(z);
  T = ncol(z);
  D_3d = D_beta(z, beta, k);
  D = drop(colSums(D_3d));
  C_3d = C_alpha(z, alpha, k);
  sum_Cj_betaj = 0;
  for (j in 1:m){
    beta_v = c(beta[j, ])
    beta_t = t(beta_v)
    beta_tt = t(beta_t)
    sum_Cj_betaj = sum_Cj_betaj + drop(C_3d[j, , ]) %*% Conj(beta_tt);
    
  }
  
  f = solve(D, sum_Cj_betaj);
  return(f);
}


F_f <- function(z, f, k){
  m = nrow(z);
  T = ncol(z);
  F = matrix(0,nrow=T,ncol=k+2)
  for (t in 1:nrow(F)){        
    F[t, ] = cbind(Conj(t(f[t:(t + k)])), 1);
  }
  return(F);
}

alpha_beta <- function(z, f, k){
  m=nrow(z)
  T=ncol(z)
  F=F_f(z, f, k)
  FtF_inv = solve(Conj(t(F)) %*% F);
  FtF_inv_Ft = FtF_inv %*% Conj(t(F));
  tmp = z %*% Conj(t(FtF_inv_Ft));
  alpha = Conj(t(tmp[, k + 2]));
  beta = tmp[, 1:(k + 1)];
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

dynamic_principal_components()

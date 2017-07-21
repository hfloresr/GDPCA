gdpc <- function(Z, k, f_ini = NULL, tol = 1e-04, niter_max = 500, crit = "LOO") {
  # A wrapper function for gdpc_priv.
  #INPUT
  # Z: data matrix each COLUMN is a different time series
  # k: number of lags used
  # f_ini: starting point for the iterations. Optional. If no argument is passed
  # the standard Principal Component completed with k leads is used.
  # tol: relative precision, stopping criterion
  # niter_max: maximum number of iterations
  # first principal component with k lags is used
  # crit: a string: "LOO", "AIC", "BIC" or "BNG"
  #OUTPUT
  # An object of class gdpc, that is, a list with entries:
  # f: coordinates of the Principal Component corresponding to the periods 1,…,T
  # initial_f: Coordinates of the Principal Component corresponding to the periods -k+1,…,0.
  # beta: beta matrix of loadings corresponding to f
  # alpha: alpha vector of intercepts corresponding to f
  # mse: mean (in T and m) squared error of the residuals of the fit
  # k: number of lags used
  # crit: the criterion of the fitted model, according to what was specified in crit
  # expart: proportion of the variance explained
  # call: the matched call
  # conv: logical. Did the iterations converge?

  sel <- switch(crit, LOO = 1, AIC = 2, BIC = 3, BNG = 4)
  if (is.null(f_ini)) {
    out <- gdpc_priv(t(Z), k, 0, FALSE, tol, niter_max, sel)
  } else {
    out <- gdpc_priv(t(Z), k, f_ini, TRUE, tol, niter_max, sel)
  }
  out$expart <- 1 - out$mse/mean(apply(Z, 2, var))
  fn_call <- match.call()
  fn_call$crit <- crit
  out$call <- fn_call
  out <- construct.gdpc(out, Z)
  return(out)
}

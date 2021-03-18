
## covarIW_helper.R

## functions included:
## (1) logmultigamma() : log multivariate gamma function
## (2) sampleIW()      : sample from IW-distribution
## (3) preprocess()    : evaluate psi(u) for each of the posterior samples
## (4) lil()           : compute true log integrated (marginal) likelihood
## (5) maxLogLik()     : compute maximized likelihood using true Sigma
## (6) cov_loglik()    : compute log likelihood
## (7) cov_logprior()  : compute log prior (IW)
## (8) psi()           : compute -loglikelihood-logprior

logmultigamma = function(p, a) {
  f = 0.25 * p * (p - 1) * log(pi)
  for(i in 1:p) { f = f + lgamma( a + 0.5 - 0.5 * i) }
  return(f)
} # end logmultigamma() function -----------------------------------------------





## getDiagIndex() function
# return the indices of the diagonal elements in the cholesky factor
# D   : dimension of the covariance matrix
# D_u : number of elements in the lower cholesky factor (dimension of u)
getDiagIndex = function(D, D_u) {

  L_ld = matrix(0, D, D)
  L_ld[lower.tri(L_ld, diag = T)] = 1:D_u
  diag_ind = diag(L_ld)

  return(diag_ind)

} # end getDiagIndex() function ------------------------------------------------





getDiagCols = function(hml_obj, D, D_u) {

  ind = getDiagIndex(D, D_u)
  diag_names = paste("u", ind, sep = '')

  dplyr::select(hml_obj, matches(diag_names))

}





## sampleIW() function
# input:
#        J     : # of MCMC samples to draw from inverse wishart
#        N     : sample size
#        D_u   : dim of vector u (lower triagular elements of cholesky factor)
#        nu    : prior degrees of freedom
#        S     : sum_n x_n x_n'
#        Omega : prior scale matrix
sampleIW = function(J, N, D_u, nu, S, Omega) {


  Sigma_post = vector("list", J)  # store posterior samples in matrix form
  L_post     = vector("list", J)  # store lower cholesky factor in matrix form
  post_samps = matrix(0, J, D_u)  # store posterior samples in vector form

  for (j in 1:J) {

    # Sigma_post_j = solve(Wishart_InvA_RNG(nu + N, S + Omega))

    Sigma_post_j = MCMCpack::riwish(nu + N, Omega + S)

    L_post_j = t(chol(Sigma_post_j))

    Sigma_post[[j]] = Sigma_post_j                             # (D x D)
    L_post[[j]]     = L_post_j                                 # (D x D)
    post_samps[j,]  = L_post_j[lower.tri(L_post_j, diag = T)]  # (D_u x 1)

  } # end of sampling loop



  return(list(post_samps = data.frame(post_samps), Sigma_post = Sigma_post,
              L_post = L_post))


} # end sampleIW() function ----------------------------------------------------





lil = function(param_list) {

  N     = param_list$N      # number of observations
  nu    = param_list$nu     # prior degrees of freedom for inv-wishart
  D     = param_list$D      # num cols/rows in covariance matrix Sigma
  Omega = param_list$Omega  # prior scale matrix for inv-wishart
  S     = param_list$S      # sum_n x_n x_n'

  logML = logmultigamma(D, (nu + N) / 2) - 0.5 * N * D * log(pi) -
    logmultigamma(D, nu / 2) + 0.5 * nu * log_det(Omega) -
    0.5 * (N + nu) * log_det(Omega + S)

  return(logML)

} # end logML() function -------------------------------------------------------





## maxLogLik() function
maxLogLik = function(Sigma, param_list) {

  N     = param_list$N      # number of observations
  D     = param_list$D      # num cols/rows in covariance matrix Sigma
  S     = param_list$S      # sum_n x_n x_n'

  loglik = - 0.5 * N * D * log(2 * pi) - 0.5 * N * log_det(Sigma) -
    0.5 * matrixcalc::matrix.trace(solve(Sigma) %*% S)

  return(loglik)

} # end maxLogLik() function ---------------------------------------------------





# cov_loglik() function -------------------------------------------------------
cov_loglik = function(u, params) {

  N = params$N
  D = params$D
  S = params$S

  L = matrix(0, D, D)             # (D x D) lower triangular matrix
  L[lower.tri(L, diag = T)] = u   # populate lower triangular terms

  logDiagL = log(diag(L))         # log of diagonal terms of L

  # compute loglikelihood using the 'L' instead of 'Sigma' --> allows us
  # to use simplified version of the determinant
  loglik = - 0.5 * N * D * log(2 * pi) - N * sum(logDiagL) -
    0.5 * matrixcalc::matrix.trace(solve(L %*% t(L)) %*% S)

  return(loglik)

} # end of cov_loglik() function ---------------------------------------------





## cov_logprior() function  --------------------------------------------------
cov_logprior = function(u, params) {

  Omega = params$Omega # prior scale matrix
  nu = params$nu       # prior degrees of freedom
  D  = params$D        # dimension of covariance matrix Sigma, L

  L = matrix(0, D, D)             # (D x D) lower triangular matrix
  L[lower.tri(L, diag = T)] = u   # populate lower triangular terms

  logDiagL = log(diag(L))         # log of diagonal terms of L

  # compute log of constant term
  logC = 0.5 * nu * log_det(Omega) - 0.5 * (nu * D) * log(2) -
    logmultigamma(D, nu / 2)

  # compute log of the Jacobian term
  logJacTerm = D * log(2) + sum((D + 1 - 1:D) * logDiagL)

  # compute full log prior expression -- note log(det(L)) = sum(logDiagL)
  logprior = logC - (nu + D + 1) * sum(logDiagL) -
    0.5 * matrixcalc::matrix.trace(solve(L %*% t(L)) %*% Omega) +
    logJacTerm

  return(logprior)

} # end of cov_logprior() function -------------------------------------------





## psi() function  -----------------------------------------------------------
psi_covar = function(u, params) {
  loglik = cov_loglik(u, params)
  logprior = cov_logprior(u, params)
  return(- loglik - logprior)
} # end of psi() function ----------------------------------------------------


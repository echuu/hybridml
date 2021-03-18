

source("examples/iw/iw_helper.R")
Rcpp::sourceCpp("examples/iw/fast_iw.cpp")

N = 100                     # number of observations
D = 6                       # num rows/cols in the covariance matrix
D_u = 0.5 * D * (D + 1)     # dimension of u that is fed into the tree
J = 1000


## wishart prior parameters
Omega = diag(1, D)          # scale matrix
nu    = D + 1               # degrees of freedom

## specify the true covariance matrix
set.seed(1)
Sigma = matrix(stats::rWishart(1, D, Omega), D)
is.positive.definite(Sigma)

## (1) generate data
X = mvtnorm::rmvnorm(N, mean = rep(0, D), sigma = Sigma) # (N x p)
S = t(X) %*% X                                  # (p x p)

## store parameters in a list that can be passed into the algorithm
params = list(S = S, N = N, D = D, D_u = D_u, # S, dimension vars
              Omega = Omega, nu = nu)         # prior params

## (2) obtain posterior samples
postIW = sampleIW(J, N, D_u, nu, S, Omega)     # post_samps, Sigma_post, L_post

# these are the posterior samples stored as vectors (the lower cholesky factors
# have been collapsed into D_u dimensional vectors)
post_samps = postIW$post_samps                 # (J x D_u)
u_df = hybridml::preprocess(post_samps, D_u, params)

## true log marginal likelihood
(LIL = lil(param_list))

lambda = function(u, params) { pracma::grad(psi_covar, u, params = params) }
hess   = function(u, params) { pracma::hessian(psi_covar, u, params = params) }

## hybrid approximation (constant)
out1 = hybridml::hybml_const(u_df)
out1$logz
out1$bounds

## hybrid approximation (EP)
out2 = hybridml::hybml(u_df, param_list, grad = lambda, hess = hess)
out2$logz
out2$bounds





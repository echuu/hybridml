
### updated gwishart example
library(BDgraph)
library(dplyr)

I_G = function(delta) {
  7/2 * log(pi) + lgamma(delta + 5/2) - lgamma(delta + 3) +
    lgamma(delta + 1) + lgamma(delta + 3/2) + 2 * lgamma(delta + 2) +
    lgamma(delta + 5/2)
}


p = 5
G_5 = matrix(c(1,1,0,1,1,
               1,1,1,0,0,
               0,1,1,1,1,
               1,0,1,1,1,
               1,0,1,1,1), p, p)
b = 3
n = 10
V = n * diag(1, p)



edgeInd = G_5[upper.tri(G_5, diag = TRUE)] %>% as.logical

# (LIL = log(2^(0.5*p*b + 7)) + I_G(0.5*(b-2)))

(truth = log(2^(0.5*p*b + 7)) + I_G(0.5*(b-2)) + (-0.5 * p * b - 7) * log(n))
gnorm(G_5, b, V, 100)

## compute some values used in the log posterior formula
nu_i = numeric(p)
for (i in 1:p) {
  ne = which(G_5[i,] > 0)
  nu_i[i] = length(ne[ne > i])
}

FREE_PARAMS_ALL = c(upper.tri(diag(1, p), diag = T) & G_5)


## construct A matrix so that we can compute k_i
A = (upper.tri(diag(1, p), diag = F) & G_5) + 0

k_i  = colSums(A) # see step 2, p. 329 of Atay
nu_i = rowSums(A) # see step 2, p. 329 of Atay
b_i = nu_i + k_i + 1
b_i

set.seed(1)
Omega_G = rgwish(1, G_5, b, V) # generate the true precision matrix

P = chol(solve(V)) # upper cholesky factor; D^(-1) = TT'  in Atay paper

params = list(G_5 = G_5, P = P, p = p, edgeInd = edgeInd,
              b = b, nu_i = nu_i, b_i = b_i)

J = 5000
N = 0
S = matrix(0, p, p)
D = sum(edgeInd) # number of free parameters / dimension of parameter space

########################

# Omega_post = rgwish(5, G_5, b + N, V + S)
#
# ## Compute Phi (upper triangular), stored column-wise, Phi'Phi = Omega
# Phi = apply(Omega_post, 3, chol) # (p^2) x J
#
# ## Compute Psi
# Psi = apply(Phi, 2, computePsi, P = P)
#
# matrix(Psi[,1], p, p)


########################

samps = samplegw(J, G_5, b, N, V, S, solve(P), FREE_PARAMS_ALL)
u_samps = samps$Psi_free %>% data.frame

u_df = preprocess(u_samps, D, params)     # J x (D_u + 1)
u_df %>% head

grad = function(u, params) { pracma::grad(psi, u, params = params) }
hess = function(u, params) { pracma::hessian(psi, u, params = params) }
u_star = globalMode(u_df)
u_star

### (1) find global mean
# MAP_LOC = which(u_df$psi_u == min(u_df$psi_u))
# u_0 = u_df[MAP_LOC,1:D] %>% unname() %>% unlist()
# u_star = u_0
gnorm(G_5, b, V, 100)
logzhat = hybridml::hybml(u_df, params, psi = psi, grad = grad, hess = hess, u_0 = u_star)$logz
logzhat
# hybridml::hybml_const(u_df)$logz

(truth - logzhat)


n_sims = 100
hyb = numeric(n_sims)
gnorm_approx = numeric(n_sims)
bridge = numeric(n_sims)
j = 1
set.seed(1)
while (j <= n_sims) {

  samps = samplegw(J, G_5, b, N, V, S, solve(P), FREE_PARAMS_ALL)
  u_samps = samps$Psi_free %>% data.frame

  u_df = preprocess(u_samps, D, params)     # J x (D_u + 1)
  # u_df %>% head

  hyb[j] = hybridml::hybml(u_df, params, psi = psi, grad = grad,
                           hess = hess, u_0 = u_star)$logz

  ### bridge estimator ---------------------------------------------------------
  u_samp = as.matrix(u_samps)
  colnames(u_samp) = names(u_df)[1:D]
  # prepare bridge_sampler input()
  lb = rep(-Inf, D)
  ub = rep(Inf, D)
  names(lb) <- names(ub) <- colnames(u_samp)

  bridge_result = bridgesampling::bridge_sampler(samples = u_samp,
                                                 log_posterior = log_density,
                                                 data = params,
                                                 lb = lb, ub = ub,
                                                 silent = TRUE)
  bridge[j] = bridge_result$logml
  ### bridge estimator ---------------------------------------------------------

  gnorm_approx[j] = gnorm(G_5, b, V, J)

  print(paste('iter ', j, ': ',
              round(mean(hyb[hyb!=0]), 3),
              ' (error = ',
              round(mean(truth - hyb[hyb!=0]), 3), ')',
              sep = ''))

  print(paste('iter ', j, ': ',
              round(mean(gnorm_approx[gnorm_approx!=0]), 3),
              ' (error = ',
              round(mean(truth - gnorm_approx[gnorm_approx!=0]), 3),
              ')', sep = ''))
  j = j + 1
}



approx = data.frame(truth, hyb = hyb, gnorm = gnorm_approx, bridge = bridge)
approx_long = reshape2::melt(approx, id.vars = 'truth')


data.frame(logz = colMeans(approx), approx_sd = apply(approx, 2, sd),
           avg_error = colMeans(truth - approx))



# test hybrid algorithm --------------------------------------------------------

## (2) fit the regression tree via rpart()
u_rpart = rpart::rpart(psi_u ~ ., u_df)

## (3) process the fitted tree
# (3.1) obtain the (data-defined) support for each of the parameters
param_support = extractSupport(u_df, D) #

# (3.2) obtain the partition
u_partition = extractPartition(u_rpart, param_support)

#### extension starts here -------------------------------------------------

log_density = function(u, data) {
  -psi(u, data)
}

### (1) find global mean
grad = function(u, params) { pracma::grad(psi, u, params = params) }
hess = function(u, params) { pracma::hessian(psi, u, params = params) }
u_star = globalMode(u_df)
u_star

### (1) find global mean
MAP_LOC = which(u_df$psi_u == min(u_df$psi_u))
u_0 = u_df[MAP_LOC,1:D] %>% unname() %>% unlist()
u_star = u_0


### (2) find point in each partition closest to global mean (for now)
# u_k for each partition
u_df_part = u_df %>% dplyr::mutate(leaf_id = u_rpart$where)

cost = apply(u_df_part[,1:D], 1, l1_norm, u_0 = u_star)
u_df_part = u_df_part %>% dplyr::mutate(cost = cost)

# take min result, group_by() leaf_id
psi_df = u_df_part %>%
  group_by(leaf_id) %>% filter(cost == min(cost)) %>%
  data.frame
# psi_df
# ------------------------------------------------------------------------------
bounds = u_partition %>% arrange(leaf_id) %>%
  dplyr::select(-c("psi_hat", "leaf_id"))
psi_df = psi_df %>% arrange(leaf_id)
# psi_df

K = nrow(bounds)
log_terms = numeric(K) # store terms so that we can use log-sum-exp()
G_k = rep(NA, K)       # store terms coming from gaussian integral


k = 1
# source("C:/Users/ericc/rcpp-epmgp/util.R")
for (k in 1:nrow(bounds)) {

  u_k = unname(unlist(psi_df[k,1:D]))

  # H_k = pracma::hessian(slow_psi, u_k, params = params)
  H_k = hess(u_k, params)
  H_k_inv = chol2inv(chol(H_k))

  # lambda_k = pracma::grad(slow_psi, u_k, params = params)
  lambda_k = grad(u_k, params = params)
  b_k = H_k %*% u_k - lambda_k
  m_k = H_k_inv %*% b_k

  lb = bounds[k, seq(1, 2 * D, 2)] %>% unname %>% unlist
  ub = bounds[k, seq(2, 2 * D, 2)] %>% unname %>% unlist

  # G_k[k] = hybridml::epmgp_stable(m_k, H_k_inv, b_k, lb, ub, EPS_CONVERGE = 1e-5)$logZ
  G_k[k] = epmgp::pmvn(lb, ub, m_k, H_k_inv, log = TRUE)
  # G_k[k] = log(TruncatedNormal::pmvnorm(m_k, H_k_inv, lb, ub)[1])

  log_terms[k] = D / 2 * log(2 * pi) - 0.5 * log_det(H_k) -
    psi_df$psi_u[k] + sum(lambda_k * u_k) - 0.5 * t(u_k) %*% H_k %*% u_k +
    0.5 * t(m_k) %*% H_k %*% m_k + G_k[k]

}

log_terms

log_sum_exp(log_terms)
gnorm(G_5, b, V, 2000)
LIL

abs(LIL - gnorm(G_5, b, V, 1000))
abs(LIL - log_sum_exp(log_terms))


#### bridge sampler
log_density = function(u, data) {
  -psi(u, data)
}

u_samp = as.matrix(u_samps)
colnames(u_samp) = names(u_df)[1:D]
# prepare bridge_sampler input()
lb = rep(-Inf, D)
ub = rep(Inf, D)
names(lb) <- names(ub) <- colnames(u_samp)

bridge_result = bridgesampling::bridge_sampler(samples = u_samp,
                                               log_posterior = log_density,
                                               data = params,
                                               lb = lb, ub = ub,
                                               silent = TRUE)
bridge_result$logml
abs(LIL - bridge_result$logml)













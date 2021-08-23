
## 5/10:
## Split the calculation in to blocks, since each block of G_5 gives a very
## easy calculation for both the gradient and the hessian.
## This will require separate functions for both the gradient and the hessian
## that are specific to this problem. In addition, we have to modify
## globalMode() so that the correct params_G5 is passed into the grad and hess
## functions, as well as hybridml(), so that the correct params_G5 is passed
## into the grad and hess functions


source("examples/gwish/gwish_density.R")
library(BDgraph)
library(dplyr)
#### initialize graphs ---------------------------------------------------------



### create parameters for a single G_5 -----------------------------------------
p1 = 5
G_5 = matrix(c(1,1,0,1,1,
               1,1,1,0,0,
               0,1,1,1,1,
               1,0,1,1,1,
               1,0,1,1,1), p1, p1)

G_5 = matrix(1,5,5)


p1 = 10
n = 3
G_5 = matrix(c(0, 1, 0, 0, 0, 0, 0, 1, 0, 0,
               1, 0, 1, 0, 0, 0, 1, 0, 0, 0,
               0, 1, 0, 1, 0, 1, 0, 0, 0, 0,
               0, 0, 1, 0, 1, 0, 0, 0, 0, 0,
               0, 0, 0, 1, 0, 1, 0, 0, 0, 0,
               0, 0, 1, 0, 1, 0, 1, 0, 0, 0,
               0, 1, 0, 0, 0, 1, 0, 1, 0, 0,
               1, 0, 0, 0, 0, 0, 1, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 10, 10)
diag(G_5) = 1



p1 = 6
G_5 = matrix(1,6,6)
G_5[6,] = 0
G_5[,6] = 0
G_5[6,6] = 1
gnorm(G_5, 3, 3*diag(6), 100)



b = 5
n = 3
V_5 = n * diag(1, p1)

P = chol(solve(V_5)) # upper cholesky factor; D^(-1) = TT'  in Atay paper

FREE_PARAMS_ALL = c(upper.tri(diag(1, p1), diag = T) & G_5)
edgeInd = G_5[upper.tri(G_5, diag = TRUE)] %>% as.logical

## construct A matrix so that we can compute k_i
A = (upper.tri(diag(1, p1), diag = F) & G_5) + 0

k_i  = colSums(A) # see step 2, p. 329 of Atay
nu_i = rowSums(A) # see step 2, p. 329 of Atay
b_i = nu_i + k_i + 1
b_i

# S = matrix(0, p, p)
D = sum(edgeInd) # number of free parameters / dimension of parameter space

index_mat = matrix(0, p1, p1)
index_mat[upper.tri(index_mat, diag = T)][edgeInd] = 1:D
index_mat[upper.tri(index_mat, diag = T)]
t_ind = which(index_mat!=0,arr.ind = T)
t_ind

index_mat[lower.tri(index_mat)] = NA
vbar = which(index_mat==0,arr.ind = T) # non-free elements
n_nonfree = nrow(vbar)

params_G5 = list(G = G_5, P = P, p = p1, D = D, edgeInd = edgeInd,
                 b = b, nu_i = nu_i, b_i = b_i,
                 t_ind = t_ind, n_nonfree = n_nonfree, vbar = vbar,
                 n_graphs = 1)


################################################################################

### create parameters for stacked G_5's so we can sample, and compute psi ------

n_G5 = 1 # number of G_5 graphs we want to stack
G = diag(1, n_G5) %x% G_5
p = ncol(G)
V = n * diag(1, p)

## try computing the normalizing constant of G_9 first as sanity check
FREE_PARAMS_ALL = c(upper.tri(diag(1, p), diag = T) & G)

edgeInd = G[upper.tri(G, diag = TRUE)] %>% as.logical

## construct A matrix so that we can compute k_i
A = (upper.tri(diag(1, p), diag = F) & G) + 0

k_i  = colSums(A) # see step 2, p. 329 of Atay
nu_i = rowSums(A) # see step 2, p. 329 of Atay
b_i = nu_i + k_i + 1
b_i

set.seed(1)
Omega_G = rgwish(1, G, b, V) # generate the true precision matrix
P = chol(solve(V)) # upper cholesky factor; D^(-1) = TT'  in Atay paper

# params = list(G = G, P = P, p = p, edgeInd = edgeInd,
#               b = b, nu_i = nu_i, b_i = b_i)
# N = 0
S = matrix(0, p, p)
D = sum(edgeInd) # number of free parameters / dimension of parameter space

index_mat = matrix(0, p, p)
index_mat[upper.tri(index_mat, diag = T)][edgeInd] = 1:D
index_mat[upper.tri(index_mat, diag = T)]
t_ind = which(index_mat!=0,arr.ind = T)
# t_ind

index_mat[lower.tri(index_mat)] = NA
vbar = which(index_mat==0,arr.ind = T) # non-free elements
n_nonfree = nrow(vbar)

params = list(G = G, P = P, p = p, D = D, edgeInd = edgeInd,
              b = b, nu_i = nu_i, b_i = b_i,
              t_ind = t_ind, n_nonfree = n_nonfree, vbar = vbar)

J = 200
samps = samplegw(J, G, b, 0, V, S, solve(P), FREE_PARAMS_ALL)
u_samps = samps$Psi_free %>% data.frame
# u_df = preprocess(u_samps, D, params)     # J x (D_u + 1)
u_df = gwish_preprocess(u_samps, D, params_G5)     # J x (D_u + 1)
u_df %>% head

I_G = function(delta) {
  7/2 * log(pi) + lgamma(delta + 5/2) - lgamma(delta + 3) +
    lgamma(delta + 1) + lgamma(delta + 3/2) + 2 * lgamma(delta + 2) +
    lgamma(delta + 5/2)
}
(Z_5  = log(2^(0.5*p1*b + 7)) + I_G(0.5*(b-2)) + (-0.5 * p1 * b - 7) * log(n))
Z = n_G5 * Z_5
Z
gnorm(G, b, V, J) # gnorm estimate of the entire (appended graph)


# ------------------------------------------------------------------------------

#### compute the hybrid approximation, bridge sampling estimator


## numerical versions of the grad and hess -- included here just as a comparison
## for smaller dimension to make sure the analytical gradient and hessian
## are computed correctly
# grad = function(u, params) { pracma::grad(psi, u, params = params) }
# hess = function(u, params) { pracma::hessian(psi, u, params = params) }
# u_star = hybridml::globalMode(u_df) ## slow version
# u_star_numer = u_star

grad = function(u, params) { fast_grad(u, params) }
hess = function(u, params) { fast_hess(u, params) }
u_star = gwish_globalMode(u_df, params, params_G5)

## laplace estimator
0.5*(D)*log(2*pi) -
  0.5*log_det(hess(u_star, params_G5)$H) -
  gwish_psi(c(u_star), params_G5)

# u_star
# u_star_closed = u_star

# cbind(u_star_numer, u_star_closed)
logzhat = hybml_gwish(u_df, params_G5, psi = psi, grad = grad, hess = hess, u_0 = u_star)$logz
logzhat = hybml_gwish(u_df, params_G5, psi = psi, grad = grad, hess = hess)$logz

logzhat           # hybrid
Z                 # truth


J = 2000
samps = samplegw(J, G, b, 0, V, S, solve(P), FREE_PARAMS_ALL)
u_samps = samps$Psi_free %>% data.frame
u_df = gwish_preprocess(u_samps, D, params_G5)     # J x (D_u + 1)
logzhat = hybml_gwish(u_df, params_G5, psi = gwish_psi, grad = grad, hess = hess, u_0 = u_star)$logz
logzhat           # hybrid

abs(logzhat - Z)


## bridge calculation
log_density = function(u, data) {
  -psi(u, data)
}
log_density = function(u, data) {
  -gwish_psi(u, data)
}


u_samp = as.matrix(u_samps)
colnames(u_samp) = names(u_df)[1:D]
# prepare bridge_sampler input()
lb = rep(-Inf, D)
ub = rep(Inf, D)
names(lb) <- names(ub) <- colnames(u_samp)

bridge_result = bridgesampling::bridge_sampler(samples = u_samp,
                                               log_posterior = log_density,
                                               data = params_G5,
                                               lb = lb, ub = ub,
                                               method = 'normal',
                                               silent = TRUE)

bridge_result = bridgesampling::bridge_sampler(samples = u_samp,
                                               log_posterior = log_density,
                                               data = params,
                                               lb = lb, ub = ub,
                                               method = 'normal',
                                               silent = TRUE)
bridge_result$logml

abs(bridge_result$logml - Z)

gnorm(G, b, V, J) # gnorm estimate of the entire (appended graph)
abs(gnorm(G, b, V, J) - Z)
abs(bridge_result$logml - Z)
abs(logzhat - Z)



### run simulations






# initialize storage for each of the estimators
n_sims       = 100
hyb          = numeric(n_sims)
# hyb_old      = numeric(n_sims)
gnorm_approx = numeric(n_sims)
bridge       = numeric(n_sims)
bridge_warp  = numeric(n_sims)

j = 1
J = 1000
set.seed(1)
while (j <= n_sims) {

  samps = samplegw(J, G, b, N, V, S, solve(P), FREE_PARAMS_ALL)
  u_samps = samps$Psi_free %>% data.frame

  ### hyb estimator ------------------------------------------------------------
  u_df = gwish_preprocess(u_samps, D, params_G5)     # J x (D_u + 1)
  logzhat = hybml_gwish(u_df, params_G5, psi = psi, grad = grad, hess = hess,
                        u_0 = u_star)$logz
  hyb[j] = logzhat # hybrid

  ### bridge estimator ---------------------------------------------------------
  u_samp = as.matrix(u_samps)
  colnames(u_samp) = names(u_df)[1:D]
  lb = rep(-Inf, D)
  ub = rep(Inf, D)
  names(lb) <- names(ub) <- colnames(u_samp)
  bridge_result = bridgesampling::bridge_sampler(samples = u_samp,
                                                 log_posterior = log_density,
                                                 data = params_G5,
                                                 lb = lb, ub = ub,
                                                 silent = TRUE)
  bridge[j] = bridge_result$logml

  bridge_result = bridgesampling::bridge_sampler(samples = u_samp,
                                                 log_posterior = log_density,
                                                 data = params_G5,
                                                 lb = lb, ub = ub,
                                                 method = 'warp3',
                                                 silent = TRUE)
  bridge_warp[j] = bridge_result$logml

  ### gnorm estimator ----------------------------------------------------------
  gnorm_approx[j] = gnorm(G, b, V, J)


  ### display some information about the running avg of the estimator + error
  print(paste('iter ', j, ': ',
              'hyb = ',     round(mean(hyb[hyb!=0]), 3),
              ' (error = ', round(mean((Z - hyb[hyb!=0])), 3), '), ',
              'bse = ',  round(mean(bridge[bridge!=0]), 3),
              ' (error = ', round(mean((Z - bridge[bridge!=0])), 3), '), ',
              'wbse = ',  round(mean((bridge_warp[bridge_warp!=0])), 3),
              ' (error = ', round(mean(Z - bridge_warp[bridge_warp!=0]), 3), '), ',
              sep = ''))

  j = j + 1
}


mean(abs(Z - hyb[hyb != 0]))
mean(abs(Z - bridge[bridge != 0]))
mean(abs(Z - bridge_warp[bridge_warp != 0]))

mean(abs(Z - hyb[hyb != 0])^2)
mean(abs(Z - bridge[bridge != 0])^2)
mean(abs(Z - bridge_warp[bridge_warp != 0])^2)

truth = Z
approx = data.frame(truth, hyb = hyb, gnorm = gnorm_approx, bridge = bridge,
                    bridge_warp = bridge_warp)
approx_long = reshape2::melt(approx, id.vars = 'truth')

truth = Z
res_tbl =
  data.frame(logz      = colMeans(approx),
             approx_sd = apply(approx, 2, sd),
             avg_error = colMeans(truth - approx),            # avg error
             mae       = colMeans(abs(truth - approx)),       # mean absolute error
             rmse      = sqrt(colMeans((truth - approx)^2)))  # root mean square error



model_name = paste("gw_", 'p_', p, '_D_', D, '_stack_', n_G5, sep = '')
results_obj = list(name = model_name,
                   p = p,
                   D = D,
                   truth = truth,
                   results = res_tbl,
                   nMCMC = J,
                   n = n,
                   delta = b,
                   all_approx = approx,
                   Omega_G = Omega_G)
save_loc = 'C:/Users/ericc/Documents/sim_results/'
paste(save_loc, model_name, ".rds", sep = '')
saveRDS(results_obj, file = paste(save_loc, model_name, ".rds", sep = ''))









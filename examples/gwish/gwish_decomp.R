
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
#### initialize graphs ---------------------------------------------------------



### create parameters for a single G_5 -----------------------------------------
p1 = 5
G_5 = matrix(c(1,1,0,1,1,
               1,1,1,0,0,
               0,1,1,1,1,
               1,0,1,1,1,
               1,0,1,1,1), p1, p1)

b = 100
n = 100
V_5 = n * diag(1, p1)

FREE_PARAMS_ALL = c(upper.tri(diag(1, p), diag = T) & G_5)
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
                 t_ind = t_ind, n_nonfree = n_nonfree, vbar = vbar)


################################################################################

### create parameters for stacked G_5's so we can sample, and compute psi ------

n_G5 = 11 # number of G_5 graphs we want to stack
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
N = 0
S = matrix(0, p, p)
D = sum(edgeInd) # number of free parameters / dimension of parameter space

index_mat = matrix(0, p, p)
index_mat[upper.tri(index_mat, diag = T)][edgeInd] = 1:D
index_mat[upper.tri(index_mat, diag = T)]
t_ind = which(index_mat!=0,arr.ind = T)
t_ind

index_mat[lower.tri(index_mat)] = NA
vbar = which(index_mat==0,arr.ind = T) # non-free elements
n_nonfree = nrow(vbar)

params = list(G = G, P = P, p = p, D = D, edgeInd = edgeInd,
              b = b, nu_i = nu_i, b_i = b_i,
              t_ind = t_ind, n_nonfree = n_nonfree, vbar = vbar)


samps = samplegw(J, G, b, N, V, S, solve(P), FREE_PARAMS_ALL)
u_samps = samps$Psi_free %>% data.frame
# u_df = preprocess(u_samps, D, params)     # J x (D_u + 1)
# u_df %>% head

I_G = function(delta) {
  7/2 * log(pi) + lgamma(delta + 5/2) - lgamma(delta + 3) +
    lgamma(delta + 1) + lgamma(delta + 3/2) + 2 * lgamma(delta + 2) +
    lgamma(delta + 5/2)
}
(Z_5  = log(2^(0.5*p1*b + 7)) + I_G(0.5*(b-2)) + (-0.5 * p1 * b - 7) * log(n))
Z = n_G5 * Z_5
Z

# ------------------------------------------------------------------------------


grad = function(u, params) { pracma::grad(psi, u, params = params) }
hess = function(u, params) { pracma::hessian(psi, u, params = params) }
u_star = hybridml::globalMode(u_df) ## slow version
u_star_numer = u_star

grad = function(u, params) { fast_grad(u, params)  }
hess = function(u, params) { fast_hess(u, params) }
u_star = gwish_globalMode(u_df, params, params_G5)
# u_star_closed = u_star

# cbind(u_star_numer, u_star_closed)
logzhat = hybml_gwish(u_df, params_G5, psi = psi, grad = grad, hess = hess, u_0 = u_star)$logz
logzhat           # hybrid
Z                 # truth


J = 1000
samps = samplegw(J, G, b, N, V, S, solve(P), FREE_PARAMS_ALL)
u_samps = samps$Psi_free %>% data.frame
u_df = preprocess(u_samps, D, params)     # J x (D_u + 1)
logzhat = hybml_gwish(u_df, params_G5, psi = psi, grad = grad, hess = hess, u_0 = u_star)$logz
logzhat           # hybrid


## bridge calculation
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
                                               method = 'normal',
                                               silent = TRUE)
bridge_result$logml

gnorm(G, b, V, J) # gnorm estimate of the entire (appended graph)
abs(gnorm(G, b, V, J) - Z)



















### updated gwishart example


p = 5
G_5 = matrix(c(1,1,0,1,1,
               1,1,1,0,0,
               0,1,1,1,1,
               1,0,1,1,1,
               1,0,1,1,1), p, p)
b = 3
V = diag(1, p)
gnorm(G_5, b, V, 100)

## compute some values used in the log posterior formula
nu_i = numeric(p)
for (i in 1:p) {
  ne = which(G_5[i,] > 0)
  nu_i[i] = length(ne[ne > i])
}

k_i = (1:p) - 1
b_i = nu_i + k_i + 1
b_i

set.seed(1)
Omega_G = rgwish(1, G_5, b, V) # generate the true precision matrix

P   = chol(solve(V)) # upper cholesky factor; D^(-1) = TT'  in Atay paper
Phi = chol(Omega_G)  # upper cholesky factor; K = Phi'Phi   in Atay paper
Psi = Phi %*% solve(P)


params = list(N = N, D = d, D_0 = D, S = S, b = b, V = V,
              G = G_5, nu = nu_i, xi = xi,
              t_ind = t_ind)

params = list(G_5 = G_5, Phi = Phi, P = P, p = p, edgeInd = edgeInd,
              b = b, nu_i = nu_i, b_i = b_i)

J = 1000
N = 0
S = matrix(0, p, p)
b = 3
D = sum(edgeInd) # number of free parameters / dimension of parameter space

samps = samplegw(J, G_5, b, N, V, S, FREE_PARAMS_ALL)
u_samps = samps$Psi_free %>% data.frame

u_df = preprocess(u_samps, D, params)     # J x (D_u + 1)
u_df %>% head

## verify that the Psi matrix formed in the log prior is the same as the Psi
## matrix with all non-free elements
u1 = u_samps[3,] %>% unname %>% unlist
psi(u1, params)
matrix(samps$Psi[3,], p, p)

cbind(pracma::grad(test, u, obj = obj),
      pracma::grad(psi,  u, params = params))


# test hybrid algorithm --------------------------------------------------------

## (2) fit the regression tree via rpart()
u_rpart = rpart::rpart(psi_u ~ ., u_df)

## (3) process the fitted tree
# (3.1) obtain the (data-defined) support for each of the parameters
param_support = extractSupport(u_df, D) #

# (3.2) obtain the partition
u_partition = extractPartition(u_rpart, param_support)

#### extension starts here -------------------------------------------------

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
psi_df
# ------------------------------------------------------------------------------
bounds = u_partition %>% arrange(leaf_id) %>%
  dplyr::select(-c("psi_hat", "leaf_id"))
psi_df = psi_df %>% arrange(leaf_id)
psi_df

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

G_k
log_terms
log_sum_exp(log_terms)
gnorm(G_5, b, V, 1000)
















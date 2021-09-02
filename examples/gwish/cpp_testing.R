

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

b = 5
n = 10
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

set.seed(1)
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




###### start hybrid method

## numerical versions of the grad and hess -- included here just as a comparison
## for smaller dimension to make sure the analytical gradient and hessian
## are computed correctly
# grad = function(u, params) { pracma::grad(psi, u, params = params) }
# hess = function(u, params) { pracma::hessian(psi, u, params = params) }
# u_star = hybridml::globalMode(u_df) ## slow version
# u_star_numer = u_star

grad = function(u, params) { fast_grad(u, params) }
hess = function(u, params) { fast_hess(u, params) }

grad = function(u, params) { grad_cpp(u, params) }
hess = function(u, params) { hess_cpp(u, params) }


u_star = gwish_globalMode(u_df, params, params_G5)

## laplace estimator
0.5*(D)*log(2*pi) -
  0.5*log_det(hess(u_star, params_G5)$H) -
  gwish_psi(c(u_star), params_G5)

# u_star
# u_star_closed = u_star

# cbind(u_star_numer, u_star_closed)
logzhat = hybml_gwish(u_df, params_G5, psi = psi, grad = grad, hess = hess, u_0 = u_star)$logz
# logzhat = hybml_gwish(u_df, params_G5, psi = psi, grad = grad, hess = hess)$logz

logzhat           # hybrid
Z                 # truth


J = 2000
samps = samplegw(J, G, b, 0, V, S, solve(P), FREE_PARAMS_ALL)
u_samps = samps$Psi_free %>% data.frame
u_df = gwish_preprocess(u_samps, D, params_G5)     # J x (D_u + 1)
logzhat = hybml_gwish(u_df, params_G5, psi = gwish_psi, grad = grad, hess = hess, u_0 = u_star)$logz

logzhat = hybml_gwish(u_df, params_G5, psi = gwish_psi, grad = grad, hess = hess, u_0 = u_star)$logz

logzhat           # hybrid


###### -------------------------------------------------------------------------

gwish_globalMode_mod = function(u_df, params, params_G5,
                            tolerance = 1e-5, maxsteps = 200) {


  # use the MAP as the starting point for the algorithm
  MAP_LOC = which(u_df$psi_u == min(u_df$psi_u))
  theta = u_df[MAP_LOC,1:params$D] %>% unname() %>% unlist()

  numsteps = 0
  tolcriterion = 100
  step.size = 1


  while(tolcriterion > tolerance && numsteps < maxsteps){
    # print(numsteps)
    # hess_obj = hess(theta, params_G5)
    G = -hess(theta, params_G5)
    invG = solve(G)
    # G = -hess(theta, params)
    # invG = solve(G)

    thetaNew = theta + step.size * invG %*% grad(theta, params_G5)

    # if precision turns negative or if the posterior probability of
    # thetaNew becomes smaller than the posterior probability of theta
    if(-psi(thetaNew, params) < -psi(theta, params)) {
      cat('tolerance reached on log scale =', tolcriterion, '\n')
      print(paste("converged -- ", numsteps, " iters", sep = ''))
      return(theta)
    }

    tolcriterion = abs(psi(thetaNew, params)-psi(theta, params))
    theta = thetaNew
    numsteps = numsteps + 1
  }

  if(numsteps == maxsteps)
    warning('Maximum number of steps reached in Newton method.')

  print(paste("converged -- ", numsteps, " iters", sep = ''))
  return(theta)
}



hybml_gwish_cpp = function(u_df, params, psi, grad, hess, u_0 = NULL, D = ncol(u_df) - 1) {

  options(scipen = 999)
  options(dplyr.summarise.inform = FALSE)

  ## fit the regression tree via rpart()
  u_rpart = rpart::rpart(psi_u ~ ., u_df)

  ## (3) process the fitted tree
  # (3.1) obtain the (data-defined) support for each of the parameters
  param_support = extractSupport(u_df, D) #

  # (3.2) obtain the partition
  u_partition = extractPartition(u_rpart, param_support)

  #### hybrid extension begins here ------------------------------------------

  ### (1) find global mean
  # u_0 = colMeans(u_df[,1:D]) %>% unname() %>% unlist() # global mean

  if (is.null(u_0)) {
    MAP_LOC = which(u_df$psi_u == min(u_df$psi_u))
    u_0 = u_df[MAP_LOC,1:D] %>% unname() %>% unlist()
    # print(u_0)
  }

  ### (2) find point in each partition closest to global mean (for now)
  # u_k for each partition
  u_df_part = u_df %>% dplyr::mutate(leaf_id = u_rpart$where)

  l1_cost = apply(u_df_part[,1:D], 1, l1_norm, u_0 = u_0)
  u_df_part = u_df_part %>% dplyr::mutate(l1_cost = l1_cost)

  # take min result, group_by() leaf_id
  psi_df = u_df_part %>%
    dplyr::group_by(leaf_id) %>% dplyr::filter(l1_cost == min(l1_cost)) %>%
    data.frame

  bounds = u_partition %>% dplyr::arrange(leaf_id) %>%
    dplyr::select(-c("psi_hat", "leaf_id"))
  psi_df = psi_df %>% dplyr::arrange(leaf_id)

  K = nrow(bounds)
  log_terms = numeric(K) # store terms so that we can use log-sum-exp()
  G_k = numeric(K)       # store terms coming from gaussian integral

  # lambda_k = apply(psi_df[,1:D], 1, lambda, params = params)
  # k = 1
  for (k in 1:K) {
    u_k = unname(unlist(psi_df[k,1:D]))

    # hess_obj = hess(u_k, params) # pass in the G5 params, NOT the G params
    H_k = hess(u_k, params)
    # H_k_inv = solve(H_k)

    # H_k = hess(u_k, params = params)
    H_k_inv = suppressWarnings(chol2inv(chol(H_k)))

    # lambda_k = pracma::grad(psi, u_k, params = params) # numerical
    lambda_k = grad(u_k, params)
    b_k = H_k %*% u_k - lambda_k
    m_k = H_k_inv %*% b_k

    lb = bounds[k, seq(1, 2 * D, 2)] %>% unname %>% unlist
    ub = bounds[k, seq(2, 2 * D, 2)] %>% unname %>% unlist

    G_k[k] = epmgp::pmvn(lb, ub, m_k, H_k_inv, log = TRUE)
    # G_k[k] = log(TruncatedNormal::pmvnorm(m_k, H_k_inv, lb, ub)[1])

    log_terms[k] = D / 2 * log(2 * pi) - 0.5 * log_det(H_k) -
      psi_df$psi_u[k] + sum(lambda_k * u_k) -
      0.5 * t(u_k) %*% H_k %*% u_k +
      0.5 * t(m_k) %*% H_k %*% m_k + G_k[k]
  }

  return(list(logz = log_sum_exp(log_terms),
              bounds = bounds,
              G_k = G_k))
}





set.seed(1)
J = 2000
samps = samplegw(J, G, b, 0, V, S, solve(P), FREE_PARAMS_ALL)
u_samps = samps$Psi_free %>% data.frame
# u_df = preprocess(u_samps, D, params)     # J x (D_u + 1)
u_df = gwish_preprocess(u_samps, D, G5_obj)     # J x (D_u + 1)
u_star = gwish_globalMode(u_df, G5_obj, G5_obj)
u_star_cpp = gwish_globalMode_mod(u_df, G5_obj, G5_obj)


hybml_gwish(u_df, params_G5, psi = psi, grad = fast_grad, hess = fast_hess, u_0 = u_star)$logz
hybml_gwish(u_df, params_G5, psi = psi, grad = fast_grad, hess = fast_hess, u_0 = u_star_cpp)$logz

library(profvis)
profvis(hybml_gwish_cpp(u_df, G5_obj, psi = psi, grad = grad_cpp, hess = hess_cpp, u_0 = u_star_cpp)$logz)
hybml_gwish_cpp(u_df, G5_obj, psi = psi, grad = grad_cpp, hess = hess_cpp, u_0 = u_star_cpp)$logz
logzhat # hybrid
Z


microbenchmark(r = hybml_gwish(u_df, params_G5, psi = psi, grad = fast_grad, hess = fast_hess, u_0 = u_star_cpp)$logz,
               cpp = hybml_gwish_cpp(u_df, G5_obj, psi = psi_cpp, grad = grad_cpp, hess = hess_cpp, u_0 = u_star_cpp)$logz,
               gnorm = gnorm(G, b, V, J),
               times = 20)


####

Rcpp::sourceCpp("C:/Users/ericc/Documents/hybridml/examples/gwish/gwish.cpp")


u = psi_df[10,1:D] %>% unlist %>% unname
log_det(hess(u, G5_obj))
approx_integral(K, as.matrix(psi_df), as.matrix(bounds), G5_obj)
approx_integral_0(K, as.matrix(psi_df), as.matrix(bounds), G5_obj)

microbenchmark(approx_integral(K, as.matrix(psi_df), G5_obj),
               approx_integral_0(K, as.matrix(psi_df), G5_obj))



u = u_df[1,1:G5_obj$D] %>% unlist %>% unname
psi_mat = create_psi_mat_cpp(u, G5_obj)

microbenchmark(r = psi(u, G5_obj),
               cpp = psi_cpp(u, G5_obj),
               cpp_faster = psi_cpp_mat(psi_mat, G5_obj))

microbenchmark(r = psi(u, G5_obj),
               cpp = psi_cpp(u, G5_obj))


hess_cpp(u, G5_obj)
ff(u, G5_obj)

microbenchmark(r = ff(u, G5_obj),
               cpp = hess_cpp(u, G5_obj))


all.equal(hess_cpp(u, G5_obj),
          ff(u, G5_obj))

all.equal(hess_cpp(u, G5_obj) %>% diag,
          ff(u, G5_obj) %>% diag)
hess_cpp(u, G5_obj) %>% diag
ff(u, G5_obj) %>% diag

ff(u, G5_obj)
ff_fast(u, G5_obj)

microbenchmark(r = ff(u, G5_obj),
               r_vec = ff_fast(u, G5_obj))


microbenchmark(r = ff(u, G5_obj),
               r_vec = ff_fast(u, G5_obj),
               cpp = grad_cpp(u, G5_obj))









#### gradient testing ----------------------------------------------------------
all.equal(create_psi_mat_cpp(u, tmp, G5_obj),
          create_psi_mat(u, G5_obj))

create_psi_mat_cpp(u, tmp, G5_obj)
create_psi_mat(u, G5_obj)
test_call(u, G5_obj)

all.equal(f(u, G5_obj),
          grad_cpp(u,  G5_obj))

microbenchmark(r = f(u, G5_obj),
               r_vec = fast_grad(u, G5_obj),
               cpp = grad_cpp(u, G5_obj))


tmp = create_psi_mat(u, G5_obj)
dpsi(1,1,tmp, G5_obj)
dpsi(1,2,tmp, G5_obj)

dpsi(1,2,tmp, G5_obj)
dpsi_cpp(0,1,tmp, G5_obj)

d1(2,4,1,2, tmp, G)
dpsi_rsij(1, 3, 0, 1, tmp, G)

u = u_df[123,1:G5_obj$D] %>% unlist %>% unname


grad_cpp(u,  G5_obj)
f(u, G5_obj)

grad_cpp(u, create_psi_mat(u, G5_obj), G5_obj)
f(u, G5_obj)
fast_grad(u, G5_obj)

library(microbenchmark)
microbenchmark(r = f(u, G5_obj),
               r_vec = fast_grad(u, G5_obj),
               cpp = grad_cpp(u, create_psi_mat_cpp(u, G5_obj), G5_obj))


test_call = function(u, params) {
  p     = params$p
  G     = params$G
  b     = params$b
  nu_i  = params$nu_i
  P     = params$P
  b_i   = params$b_i

  FREE_PARAM_MAT = upper.tri(diag(1, p), diag = T) & G
  u_mat = FREE_PARAM_MAT
  u_mat = matrix(0, p, p)
  u_mat[FREE_PARAM_MAT] = u
  u_mat
}

microbenchmark(r = test_call(u, G5_obj),
               cpp = create_psi_mat_cpp(u, tmp, G5_obj))









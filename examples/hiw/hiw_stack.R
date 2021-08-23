

G_9 = matrix(  c(1,1,0,0,1,0,0,0,0,
                 1,1,1,1,1,0,0,0,0,
                 0,1,1,1,0,0,0,0,0,
                 0,1,1,1,1,1,1,0,0,
                 1,1,0,1,1,1,0,0,0,
                 0,0,0,1,1,1,1,1,1,
                 0,0,0,1,0,1,1,1,1,
                 0,0,0,0,0,1,1,1,1,
                 0,0,0,0,0,1,1,1,1), 9, 9)

a = c(1, 3, 2, 5, 4, 6, 7, 8, 9)
G_9 = G_9[a, a]

D = nrow(G_9)
b = 3          # prior degrees of freedom
V = diag(1, D) # prior scale matrix

D_0 = 0.5 * D * (D + 1) # num entries on diagonal and upper diagonal
J = 1000

# logical vector determining existence of edges between vertices
edgeInd = testG[upper.tri(testG, diag = TRUE)] %>% as.logical
upperInd = testG[upper.tri(testG)] %>% as.logical
D_u = sum(edgeInd)


################################################################################

n_G9 = 18
G = diag(1, n_G9) %x% G_9
p = ncol(G)
V = diag(1, p)

D = nrow(G)
b = 3          # prior degrees of freedom
V = diag(1, D) # prior scale matrix

D_0 = 0.5 * D * (D + 1) # num entries on diagonal and upper diagonal
J = 1000

# logical vector determining existence of edges between vertices
edgeInd = G[upper.tri(G, diag = TRUE)] %>% as.logical
upperInd = G[upper.tri(G)] %>% as.logical
D_u = sum(edgeInd)


set.seed(1)
true_params = HIWsim(G, b, V)
Sigma_G = true_params$Sigma
Omega_G = true_params$Omega # precision matrix -- this is the one we work with


# chol(Omega_G)
# Generate data Y based on Sigma_G
N = 100
Y = matrix(0, N, D)
for (i in 1:N) {
  Y[i, ] = t(t(chol(Sigma_G)) %*% rnorm(D, 0, 1)) # (500 x D)
}

S = t(Y) %*% Y


nu = rowSums(chol(Omega_G) != 0) - 1
xi = b + nu - 1
t_ind = which(chol(Omega_G) != 0, arr.ind = T)


params = list(N = N, D = D, D_0 = D_0, D_u = D_u,
              testG = G, edgeInd = edgeInd,
              upperInd = upperInd, S = S, V = V, b = b,
              nu = nu, xi = xi, G = G, t_ind = t_ind)

################################################################################

postIW = sampleHIW(J, D_u, D_0, G, b, N, V, S, edgeInd)
post_samps = postIW$post_samps                 # (J x D_u)
u_df = hybridml::preprocess(post_samps, D_u, params)     # J x (D_u + 1)

## truth
(LIL = logmarginal(Y, G, b, V, S))

## gnorm
- 0.5 * D * N * log(2 * pi) + BDgraph::gnorm(G, b + N, V + S, iter = 5000) -
  BDgraph::gnorm(G, b, V, iter = 5000)

u_star = globalMode(u_df, params)
out = hybridml::hybml(u_df, params, grad = grad, hess = hess, u_0 = u_star)
# out = hybridml::hybml(u_df, params, grad = fast_grad, hess = fast_hess, u_0 = u_star)
out$logz


log_density = function(u, data) {
  -old_psi(u, data)
}
u_samp = as.matrix(post_samps)
colnames(u_samp) = names(u_df)[1:D_u]
# prepare bridge_sampler input()
lb = rep(-Inf, D_u)
ub = rep(Inf, D_u)
names(lb) <- names(ub) <- colnames(u_samp)

bridge_result = bridgesampling::bridge_sampler(samples = u_samp,
                                               log_posterior = log_density,
                                               data = params,
                                               lb = lb, ub = ub,
                                               silent = TRUE)
bridge_result$logml



###### run simulations ---------------------------------------------------------

postIW = sampleHIW(J, D_u, D_0, G, b, N, V, S, edgeInd)
post_samps = postIW$post_samps                 # (J x D_u)
u_df = hybridml::preprocess(post_samps, D_u, params)     # J x (D_u + 1)
u_star = globalMode(u_df, params)

# laplace
0.5*(D_u)*log(2*pi) - 0.5*log_det(hess(u_star, params)) - psi(u_star, params)
(LIL = logmarginal(Y, G, b, V, S))

n_sims       = 100
hyb          = numeric(n_sims)
bridge       = numeric(n_sims)
bridge_warp  = numeric(n_sims)
Z = LIL

j = 1
J = 1000
set.seed(1)
while (j <= n_sims) {

  postIW = sampleHIW(J, D_u, D_0, G, b, N, V, S, edgeInd)
  post_samps = postIW$post_samps                 # (J x D_u)
  u_df = hybridml::preprocess(post_samps, D_u, params)     # J x (D_u + 1)

  ### hyb estimator ----------------------------------------------------------

  logzhat = hybridml::hybml(u_df, params, grad = grad, hess = hess, u_0 = u_star)$logz
  hyb[j] = logzhat # hybrid

  ### bridge estimator -------------------------------------------------------
  u_samp = as.matrix(post_samps)
  colnames(u_samp) = names(u_df)[1:D_u]
  lb = rep(-Inf, D_u)
  ub = rep(Inf, D_u)
  names(lb) <- names(ub) <- colnames(u_samp)

  bridge_result = bridgesampling::bridge_sampler(samples = u_samp,
                                                 log_posterior = log_density,
                                                 data = params,
                                                 lb = lb, ub = ub,
                                                 silent = TRUE)
  bridge[j] = bridge_result$logml

  bridge_result = bridgesampling::bridge_sampler(samples = u_samp,
                                                 log_posterior = log_density,
                                                 data = params,
                                                 lb = lb, ub = ub,
                                                 method = 'warp3',
                                                 silent = TRUE)
  bridge_warp[j] = bridge_result$logml

  ### gnorm estimator --------------------------------------------------------
  # gnorm_approx[j] = gnorm(G, b, V, J)


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



###### compute results -------------------------------------------------------

truth = LIL
approx = data.frame(truth, hyb = hyb, bridge = bridge,
                    bridge_warp = bridge_warp)
approx = data.frame(truth = truth, hyb = hyb, bridge = bridge, bridge_warp = bridge_warp)

approx_long = reshape2::melt(approx, id.vars = 'truth')

res_tbl =
  data.frame(logz      = colMeans(approx),
             approx_sd = apply(approx, 2, sd),
             avg_error = colMeans(truth - approx),            # avg error
             mae       = colMeans(abs(truth - approx)),       # mean absolute error
             rmse      = sqrt(colMeans((truth - approx)^2)))  # root mean square error
res_tbl

t(res_tbl)[-c(2,4), c(1,3,4,2)]





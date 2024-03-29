

## create a graph object constructor to avoid too much clutter in the graphical
## model examples


init_graph = function(G, b, n, V) {

    # input:
    # G: graph in adjacency matrix form
    # b: degrees of freedom for g-wishart distribution
    # n: pseudo-sample size, used to initialize the scale matrix
    # V: unscaled matrix
    #
    # output: list with various quantities that will be used in future calcs
    # G:         graph
    # b:         degrees of freedom
    # n:         pseudo-sample size, used to initialize the scale matrix
    # p:         dimension of graph
    # V_n:       unscaled matrix
    # S:         S = X'X when there is data (for posteriors)
    # P:         upper cholesky factor of the scale matrix V_n
    # D:         dimension of parameter space (# of free parameters)
    # edgeInd:   bool indicator for edges in G (upper tri + diag elements only)
    # nu_i       defined in step 2, p. 329 of Atay
    # b_i:       defined in step 2, p. 329 of Atay
    # t_ind:     (row,col) of the location of each of the free parameters
    # n_nonfree: number of nonfree parameters
    # vbar:      nonfree elements

    p   = ncol(G)           # dimension fo the graph
    V_n = n * V             # scale matrix for gW distribution
    P   = chol(solve(V_n))  # upper cholesky factor; D^(-1) = TT' in Atay paper
    S   = matrix(0, p, p)   # S = X'X when there is data (for posteriors)

    FREE_PARAMS_ALL = c(upper.tri(diag(1, p), diag = T) & G)
    edgeInd = G[upper.tri(G, diag = TRUE)] %>% as.logical

    ## construct A matrix so that we can compute k_i
    A = (upper.tri(diag(1, p), diag = F) & G) + 0

    k_i  = colSums(A) # see step 2, p. 329 of Atay
    nu_i = rowSums(A) # see step 2, p. 329 of Atay
    b_i = nu_i + k_i + 1

    D = sum(edgeInd) # number of free parameters / dimension of parameter space

    index_mat = matrix(0, p, p)
    index_mat[upper.tri(index_mat, diag = T)][edgeInd] = 1:D
    index_mat[upper.tri(index_mat, diag = T)]
    t_ind = which(index_mat!=0,arr.ind = T)

    index_mat[lower.tri(index_mat)] = NA
    vbar = which(index_mat==0,arr.ind = T) # non-free elements
    n_nonfree = nrow(vbar)

    obj = list(G = G, b = b, n = n, p = p, V_n = V_n, S = S, P = P,
               D = D, edgeInd = edgeInd, FREE_PARAMS_ALL = FREE_PARAMS_ALL,
               nu_i = nu_i, b_i = b_i,
               t_ind = t_ind, n_nonfree = n_nonfree, vbar = vbar,
               n_graphs = 1)

    return(obj)
}


init_stacked_graph = function(G_in, b, n, V, n_graphs) {

    # input:
    # G_in:      graph in adjacency matrix form
    # b:         degrees of freedom for g-wishart distribution
    # n:         pseudo-sample size, used to initialize the scale matrix
    # V:         unscaled matrix
    # n_graphs:  number of times to stack the input graph G_in
    #
    # output:    list with various quantities that will be used in future calcs
    # G:         graph
    # b:         degrees of freedom
    # n:         pseudo-sample size, used to initialize the scale matrix
    # p:         dimension of graph
    # V_n:       unscaled matrix
    # S:         S = X'X when there is data (for posterior distributions)
    # P:         upper cholesky factor of the scale matrix V_n
    # D:         dimension of parameter space (# of free parameters)
    # edgeInd:   bool indicator for edges in G (upper tri + diag elements only)
    # nu_i       defined in step 2, p. 329 of Atay
    # b_i:       defined in step 2, p. 329 of Atay
    # t_ind:     (row,col) of the location of each of the free parameters
    # n_nonfree: number of nonfree parameters
    # vbar:      nonfree elements
    #

    G   = diag(1, n_graphs) %x% G_in ## stack the original graph n_graph times
    p   = ncol(G)                    ## dimension of the new graph
    V_n = n * V                      ## scale matrix for gW distribution
    P = chol(solve(V_n))             ## upper cholesky factor; D^(-1) = TT'
    S = matrix(0, p, p)

    FREE_PARAMS_ALL = c(upper.tri(diag(1, p), diag = T) & G)
    edgeInd = G[upper.tri(G, diag = TRUE)] %>% as.logical

    ## construct A matrix so that we can compute k_i
    A    = (upper.tri(diag(1, p), diag = F) & G) + 0
    k_i  = colSums(A) # see step 2, p. 329 of Atay
    nu_i = rowSums(A) # see step 2, p. 329 of Atay
    b_i  = nu_i + k_i + 1

    # Omega_G = rgwish(1, G, b, V) # generate the true precision matrix
    # S = matrix(0, p, p)
    D = sum(edgeInd) # number of free parameters / dimension of parameter space

    index_mat = matrix(0, p, p)
    index_mat[upper.tri(index_mat, diag = T)][edgeInd] = 1:D
    index_mat[upper.tri(index_mat, diag = T)]
    t_ind = which(index_mat!=0,arr.ind = T)

    index_mat[lower.tri(index_mat)] = NA
    vbar = which(index_mat==0,arr.ind = T) # non-free elements
    n_nonfree = nrow(vbar)

    obj = list(G = G, b = b, n = n, p = p, V_n = V_n, S = S, P = P,
               D = D, edgeInd = edgeInd, FREE_PARAMS_ALL = FREE_PARAMS_ALL,
               nu_i = nu_i, b_i = b_i,
               t_ind = t_ind, n_nonfree = n_nonfree, vbar = vbar,
               n_graphs = n_graphs)

    return(obj)

}


h = function() {
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
  approx_integral(K, as.matrix(psi_df), as.matrix(bounds), G5_obj)
}





source("examples/gwish/gwish_density.R")
library(BDgraph)
library(dplyr)

grad = function(u, params) { fast_grad(u, params) }
hess = function(u, params) { fast_hess(u, params) }

G_5 = matrix(c(1,1,0,1,1,
               1,1,1,0,0,
               0,1,1,1,1,
               1,0,1,1,1,
               1,0,1,1,1), 5, 5)

G = G_5
n = 10
V = diag(1, ncol(G_5))
b = 5
J = 200
p = ncol(G_5)

G5_obj = init_graph(G = G_5, b = 5, n = n, V = diag(1, ncol(G_5)))

I_G = function(delta) {
  7/2 * log(pi) + lgamma(delta + 5/2) - lgamma(delta + 3) +
    lgamma(delta + 1) + lgamma(delta + 3/2) + 2 * lgamma(delta + 2) +
    lgamma(delta + 5/2)
}
(Z_5  = log(2^(0.5*p*b + 7)) + I_G(0.5*(b-2)) + (-0.5 * p * b - 7) * log(n))
Z = G5_obj$n_graphs * Z_5

# G5_obj$vbar
# G5_obj$G
# G5_obj$n_graphs

set.seed(1)
J = 2000
samps = samplegw(J, G5_obj$G, G5_obj$b, 0, G5_obj$V_n, G5_obj$S, solve(G5_obj$P),
                 G5_obj$FREE_PARAMS_ALL)

u_samps = samps$Psi_free %>% data.frame
u_df = gwish_preprocess(u_samps, G5_obj$D, G5_obj)     # J x (D_u + 1)
u_star = gwish_globalMode(u_df, G5_obj, G5_obj)

logzhat = hybml_gwish(u_df, G5_obj,
                      psi = psi, grad = grad, hess = hess,
                      u_0 = u_star)$logz
logzhat # -22.72094
gnorm(G, b, n * V, J) # -22.8673
Z # -22.86793
h()

library(microbenchmark)
microbenchmark(r = hybml_gwish(u_df, params_G5, psi = psi, grad = fast_grad, hess = fast_hess, u_0 = u_star_cpp)$logz,
               cpp = hybml_gwish_cpp(u_df, G5_obj, psi = psi_cpp, grad = grad_cpp, hess = hess_cpp, u_0 = u_star_cpp)$logz,
               cpp_fast = h(),
               gnorm = gnorm(G, b, V, J),
               times = 20)








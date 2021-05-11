
#### true ML calculation
I_G = function(delta) {
  7/2 * log(pi) + lgamma(delta + 5/2) - lgamma(delta + 3) +
    lgamma(delta + 1) + lgamma(delta + 3/2) + 2 * lgamma(delta + 2) +
    lgamma(delta + 5/2)
}
(Z_5  = log(2^(0.5*p1*b + 7)) + I_G(0.5*(b-2)) + (-0.5 * p1 * b - 7) * log(n))


J = 200
gnorm(G, b, V, J) # gnorm estimate of the entire (appended graph)
(Z = n_G5 * Z_5)



grad = function(u, params) { f(u, params)  }
hess = function(u, params) { fast_hess(u, params) }
u_star = globalMode(u_df, params_G5, params)
samps = samplegw(J, G, b, N, V, S, solve(P), FREE_PARAMS_ALL)
u_samps = samps$Psi_free %>% data.frame
u_df = preprocess(u_samps, D, params)     # J x (D_u + 1)
logzhat = hybridml::hybml(u_df, params, psi = psi, grad = grad, hess = hess, u_0 = u_star)$logz
logzhat           # hybrid
Z                 # truth

####

# compute the hessian block by block

u_df %>% head
u = u_df[1,1:D] %>% unlist %>% unname
u



## TODO: also need to update hybridml that takes params_G5 as an input as well
## because params$D != params_G5$D
hybml_gwish = function(u_df, params, psi, grad, hess, u_0 = NULL, D = ncol(u_df) - 1) {

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

    hess_obj = hess(u_k, params) # pass in the G5 params, NOT the G params
    H_k = hess_obj$H
    H_k_inv = hess_obj$H_inv

    # H_k = hess(u_k, params = params)
    # H_k_inv = suppressWarnings(chol2inv(chol(H_k)))

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







## modified global mode call
gwish_globalMode = function(u_df, params, params_G5,
                            tolerance = 1e-5, maxsteps = 200) {


  # use the MAP as the starting point for the algorithm
  MAP_LOC = which(u_df$psi_u == min(u_df$psi_u))
  theta = u_df[MAP_LOC,1:params$D] %>% unname() %>% unlist()

  numsteps = 0
  tolcriterion = 100
  step.size = 1


  while(tolcriterion > tolerance && numsteps < maxsteps){
    # print(numsteps)
    hess_obj = hess(theta, params_G5)
    G = -hess_obj$H
    invG = -hess_obj$H_inv
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




fast_grad = function(u, params_G5) {

  stride = params_G5$D # dimension of parameter in the 1 graph example
  block = 1 # index for which graph we are on
  grad_list = vector("list", length = n_G5)
  for (n in 1:n_G5) {

    # extract the first (p1 x p1) block out psi_mat
    start = (block-1) * stride + 1
    end   = block * stride
    u_n = u[start:end]
    # psi_mat_n = create_psi_mat(u_n, params_G5)
    # print(psi_mat_n)
    grad_list[[n]] = f(u_n, params_G5)
    block = block + 1
  }

  return(unlist(grad_list))
}



fast_hess = function(u, params_G5) {

  stride = params_G5$D # dimension of parameter in the 1 graph example
  block = 1 # index for which graph we are on
  hess_list = vector("list", length = n_G5)
  inv_hess_list  = vector("list", length = n_G5)
  for (n in 1:n_G5) {

    # extract the first (p1 x p1) block out psi_mat
    start = (block-1) * stride + 1
    end   = block * stride
    u_n = u[start:end]
    # psi_mat_n = create_psi_mat(u_n, params_G5)
    # print(psi_mat_n)
    block = block + 1

    hess_list[[n]]     = ff_fast(u_n, params_G5)
    inv_hess_list[[n]] = solve(hess_list[[n]])
  }

  # return these
  H = matrix(Matrix::bdiag(hess_list), D, D)
  H_inv  = matrix(Matrix::bdiag(inv_hess_list), D, D)


  return(list(H = H, H_inv = H_inv))
}


create_psi_mat(u, params)
stride = params_G5$D # dimension of parameter in the 1 graph example
block = 1 # index for which graph we are on
hess_list = vector("list", length = n_G5)
inv_hess_list  = vector("list", length = n_G5)
for (n in 1:n_G5) {

  # extract the first (p1 x p1) block out psi_mat
  start = (block-1) * stride + 1
  end   = block * stride
  u_n = u[start:end]
  # psi_mat_n = create_psi_mat(u_n, params_G5)
  # print(psi_mat_n)
  block = block + 1

  hess_list[[n]]     = ff_fast(u_n, params_G5)
  inv_hess_list[[n]] = solve(hess_list[[n]])
}

# return these
H      = matrix(Matrix::bdiag(hess_list), D, D)
H_inv  = matrix(Matrix::bdiag(inv_hess_list), D, D)

big_hess      = ff_fast(u, params)
big_hess_inv  = solve(big_hess)

all.equal(H, big_hess)
all.equal(H_inv, big_hess_inv)


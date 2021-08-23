

#### 8/23:
#### TODO: documentation for this file needs to be updated, need to maybe review
####       some of these functions (analytically) and begin porting these
####       calculations to C++ for some speedup


#### 5/10: This file contains the adapted, example-specific code to calculate
#### the hybrid approximation for the gwish() model

#### uncomment next few lines to test the fast_hess(), fast_grad() function
#### on an individual data point
### u_df %>% head
### u = u_df[1,1:D] %>% unlist %>% unname
### u
### psi(u, params)


#### Functions below are specific to the gwish() problem and must be loaded
#### before attempting any sort of hybrid approximation. Also, gradient.R and
#### hessian.R contain the bulk of the implmentation for the wrapper functions
#### that are called below.


#### functions include:
#### (1) hybml_gwish()
#### (2) gwish_globalMode()
#### (3) fast_hess()
#### (4) fast_grad()

#### ---------------------------------------------------------------------------



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
  n_graphs = params_G5$n_graphs
  block = 1 # index for which graph we are on
  grad_list = vector("list", length = n_graphs)
  for (n in 1:n_graphs) {

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



fast_hess = function(u, params) {

  stride = params$D # dimension of parameter in the 1 graph example
  n_graphs = params$n_graphs
  block = 1 # index for which graph we are on
  hess_list = vector("list", length = n_graphs)
  inv_hess_list  = vector("list", length = n_graphs)
  for (n in 1:n_graphs) {

    # extract the first (p1 x p1) block out psi_mat
    start = (block-1) * stride + 1
    end   = block * stride
    u_n = u[start:end]
    # psi_mat_n = create_psi_mat(u_n, params_G5)
    # print(psi_mat_n)
    block = block + 1

    hess_list[[n]]     = ff_fast(u_n, params)
    inv_hess_list[[n]] = solve(hess_list[[n]])
  }

  # return these
  H = matrix(Matrix::bdiag(hess_list), params$D, params$D)
  H_inv  = matrix(Matrix::bdiag(inv_hess_list), params$D, params$D)


  return(list(H = H, H_inv = H_inv))
}


# psi(u, params)
# gwish_psi(u, params_G5)
#
# u_split = split(u, ceiling(seq_along(u) / 12))
# sum(unlist(lapply(u_split, psi, params = params_G5)))
#
# microbenchmark(psi       = psi(u, params),
#                psi_block = gwish_psi(u, params_G5),
#                psi_fast  = gwish_psi_fast(u, params_G5))

## pass in the smaller params
gwish_psi = function(u, params) {

  stride = params$D # dimension of parameter in the 1 graph example
  n_graphs = params$n_graphs
  block = 1 # index for which graph we are on
  psi_vec = numeric(n_graphs)

  for (n in 1:n_graphs) {
    start = (block-1) * stride + 1
    end   = block * stride
    u_n = u[start:end]
    psi_vec[n] = psi(u_n, params)
    block = block + 1
  }
  return(sum(psi_vec))
}

gwish_psi_fast = function(u, params) {
  stride = params$D
  u_split = split(u, ceiling(seq_along(u) / stride))
  sum(unlist(lapply(u_split, psi, params = params)))
}


# smaller params
gwish_preprocess = function(post_samps, D, params) {

  psi_u = apply(post_samps, 1, gwish_psi, params = params) %>% unname() # (J x 1)

  # (1.2) name columns so that values can be extracted by partition.R
  u_df_names = character(D + 1)
  for (d in 1:D) {
    u_df_names[d] = paste("u", d, sep = '')
  }
  u_df_names[D + 1] = "psi_u"

  # populate u_df
  u_df = cbind(post_samps, psi_u) # J x (D + 1)
  names(u_df) = u_df_names

  return(u_df)
}



#### ---------------------------------------------------------------------------


#### gradient related functions



#### 5/10: this is the final form of the gwish() gradient -- f()
#### load this entire file before computing any hybrid estimates
#### note: these functions are called on block matrices by the fast_grad()
#### function in gwish_calc.R

# u_df %>% head
# u = u_df[1,1:D] %>% unlist %>% unname

# psi_mat = create_psi_mat(u, params)
# psi_mat
#
# f(u, params)
# grad(u, params)
#
# f(u, params)
# grad(u, params)

f = function(u, params) {

  G = params$G
  p = params$p
  edgeInd = params$edgeInd
  psi_mat = create_psi_mat(u, params)

  # obj = list(G = G, p = p, psi_mat = psi_mat)
  gg = matrix(0, p, p)
  for (i in 1:p) {
    for (j in i:p) {
      if (G[i,j] != 0) { ## dpsi / dpsi_ij
        gg[i,j] = dpsi(i, j, psi_mat, params)
      }
    }
  }

  return(gg[upper.tri(gg, diag = TRUE)][edgeInd])
}

dpsi = function(i, j, psi_mat, params) {

  G = params$G
  p = params$p
  # psi_mat = obj$psi_mat

  if (G[i,j] == 0) { ## need this check even though it's in the main fcn b/c
    return(0)        ## we may eventually run into non-free elements
  }
  if (i == j) {

    d_ij = 0
    ### summation over non-free elements
    for (r in 1:p) {
      for (s in r:p) {
        if (G[r,s] == 0) {
          # print(paste('r = ', r, ', ', 's = ', s, ', ',
          #             'i = ', i, ', ', 'j = ', j, sep = ''))
          ## if psi_rs == 0, then no need to compute the derivative
          if (psi_mat[r, s] == 0) {
            # print("skipping derivative calculation")
            next
          }
          d_ij = d_ij + psi_mat[r,s] * d1(r, s, i, j, psi_mat, G)
        }
      }
    }

    return(d_ij - (params$b + params$nu_i[i] - 1) / psi_mat[i,i] + psi_mat[i,i])
  } else {
    d_ij = 0
    ## iterate through each entry (r,s) and compute d psi_rs / d psi_ij
    for (r in 1:p) {
      for (s in r:p) {
        if (G[r,s] == 0) {
          # print(paste('r = ', r, ', ', 's = ', s, ', ',
          #             'i = ', i, ', ', 'j = ', j, sep = ''))
          ## if psi_rs == 0, then no need to compute the derivative
          if (psi_mat[r, s] == 0) {
            # print("skipping derivative calculation")
            next
          }
          d_ij = d_ij + psi_mat[r,s] * d1(r, s, i, j, psi_mat, G)
        }
      } # end loop over s
    } # end loop over r

  } # end else

  return((d_ij + psi_mat[i,j]))
}



d1 = function(r, s, i, j, psi_mat, G) {

  # G = obj$G
  # p = obj$p
  # psi_mat = obj$psi_mat

  ## if (r,s) \in V
  if (G[r, s] > 0) {
    if (r == i && s == j) { # d psi_{ij} / d psi_{ij} = 1
      return(1)
    } else {                # d psi_{rs} / d psi_{ij} = 0, since psi_rs is free
      return(0)
    }
  }


  if (i > r)                            { return(0) } # (i,j) comes after (r,s)
  if (i == r && j > s)                  { return(0) } # same row, but after
  if (i == r && j == s && G[r, s] > 0)  { return(1) } # redundant check?
  if (i == r && j == s && G[r, s] == 0) { return(0) } # deriv wrt to non-free: 0

  if (r == 1 && s > r) { return(0) } # first row case has simplified formula

  if (r > 1) { # derivative of psi_rs in 2nd row onward

    tmp_sum = numeric(r - 1)

    for (k in 1:(r - 1)) { # deriv taken wrt Eq. (31) in Atay

      ## TODO: check values of psi_ks, psi_kr -- if 0, then can save calculation
      ##       on the derivative
      if (psi_mat[k,s] == 0 && psi_mat[k,r] == 0) {
        # print("both terms 0, save recursion")
        tmp_sum[k] = 0
        next
      } else {

        if (psi_mat[k,s] == 0) { ## if psi_ks is 0
          # print("saving on ks = 0")
          tmp_sum[k] = -1/psi_mat[r,r] *
            d1(k, s, i, j, psi_mat, G) * psi_mat[k, r]
        } else if (psi_mat[k,r] == 0) { ## if psi_kr is 0
          # print("saving on kr = 0")
          tmp_sum[k] = -1/psi_mat[r,r] *
            d1(k, r, i, j, psi_mat, G) * psi_mat[k, s]
        } else {

          if (i == j && r == i && G[r,s] == 0) {
            tmp_sum[k] = 1/psi_mat[r,r]^2 * psi_mat[k,r] * psi_mat[k,s] -
              1/psi_mat[r,r] * (
                d1(k, r, i, j, psi_mat, G) * psi_mat[k, s] +
                  d1(k, s, i, j, psi_mat, G) * psi_mat[k, r]
              )
          } else {
            tmp_sum[k] = -1/psi_mat[r,r] * (
              d1(k, r, i, j, psi_mat, G) * psi_mat[k, s] +
                d1(k, s, i, j, psi_mat, G) * psi_mat[k, r]
            )
          }
        }
      }


    } # end for()

    return(sum(tmp_sum)) ## expression derived from Eq. (99)

  } else {
    print("case not accounted for")
    return(NA)
  }

}


create_psi_mat = function(u, params) {

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

  ## u_mat should have all the free elements
  for (i in 1:p) {
    for (j in i:p) {
      if (G[i,j] > 0) {
        next # u_mat[i,j] already has value from input
      } else {
        if (i == 1) {
          u_mat[i,j] = -1/P[j,j] * sum(u_mat[i,i:(j-1)] * P[i:(j-1), j])
        } else {
          # for rows other than first row
          x0 = -1/P[j,j] * sum(u_mat[i,i:(j-1)] * P[i:(j-1), j])

          # old formulation that still depends on Phi
          # Psi_copy[i,j] = x0 -
          #   sum(Phi[1:(i-1),i] / P[i,i] / (Phi[i,i] / P[j,j]) *
          #         Phi[1:(i-1), j] / P[j,j])

          # new formulation that uses Eq. (31) in Atay -- no longer uses any
          # non-free parameters: only uses recursively defined parameters
          # and elements of the cholesky factorization of the scale matrix
          tmp = numeric(i-1)
          for (r in 1:(i-1)) {
            tmp1 = u_mat[r,i] + sum(u_mat[r,r:(i-1)] * P[r:(i-1),i] / P[i,i])
            tmp2 = u_mat[r,j] + sum(u_mat[r,r:(j-1)] * P[r:(j-1),j] / P[j,j])
            # print(length(u_mat[r,r:(j-1)]))
            # print(length(P[r:(j-1),j]))
            # print(tmp2)
            tmp[r] = tmp1 * tmp2
          }
          u_mat[i,j] = x0 - 1 / u_mat[i,i] * sum(tmp)
        }
      }
    }
  }

  return(u_mat)

}





#### hessian related functions


#### 5/10: this is the final form of the gwish hessian -- ff_fast()
#### load this entire file before computing any hybrid estimates
#### note: these functions are called on block matrices by the fast_hess()
#### function in gwish_calc.R


# index_mat = matrix(0, p, p)
# index_mat[upper.tri(index_mat, diag = T)][edgeInd] = 1:D
# index_mat[upper.tri(index_mat, diag = T)]
# t_ind = which(index_mat!=0,arr.ind = T)
# t_ind
#
# index_mat[lower.tri(index_mat)] = NA
# vbar = which(index_mat==0,arr.ind = T) # non-free elements
# n_nonfree = nrow(vbar)


# params = list(G = G_5, P = P, p = p, D = D, edgeInd = edgeInd,
#               b = b, nu_i = nu_i, b_i = b_i,
#               n_nonfree = n_nonfree, vbar = vbar, t_ind = t_ind)
#
#
# hess(u, params)
# ff(u, params)

# install.packages("profvis")
# library(profvis)
# profvis(ff(u, params))



# testblock = matrix(0, 6, 6)
# D = ncol(testblock)
# stride = 2
# block = 1
# for (r in 1:(D-1)) {
#
#   for (c in (r+1):D) {
#     if (c <= (block * stride)) {
#       testblock[r,c] = NA
#       testblock[c,r] = NA
#     } else {
#       break
#     }
#   }
#
#   if ((r %% stride) == 0) {
#     block = block + 1
#   }
#
# }
# testblock

ff_fast = function(u, params) {

  G = params$G
  D = params$D
  t_ind = params$t_ind
  n_nonfree = params$n_nonfree
  vbar = params$vbar
  psi_mat = create_psi_mat(u, params)
  #
  #
  # obj = list(G = G, p = p, psi_mat = psi_mat, n_nonfree = n_nonfree,
  #            vbar = vbar)

  H = matrix(0, D, D)
  for (d in 1:D) {

    i = t_ind[d,1] # row index
    j = t_ind[d,2] # col index

    # populate H[i,i] element
    if (i == j) { # one of the original diag elements in the psi matrix, Psi_ii
      tmp = 0
      for (a in 1:n_nonfree) {
        rr = vbar[a, 1] # row index of non-free element
        ss = vbar[a, 2] # col index of non-free element

        tmp = tmp + d1(rr, ss, i, j, psi_mat, G)^2
      }
      # if (tmp != 0) {print("tmp != 0")}

      H[d,d] = (params$b + params$nu_i[i] - 1) / psi_mat[i,i]^2 + 1 + tmp
    } else { # 2nd order partial for original off-diag elements in the psi mat

      tmp = 0
      for (a in 1:n_nonfree) {
        rr = vbar[a, 1] # row index of non-free element
        ss = vbar[a, 2] # col index of non-free element

        tmp = tmp + d1(rr, ss, i, j, psi_mat, G)^2
      }
      H[d,d] = 1 + tmp
    }
  }


  ### Populate the hessian matrix. This is a (D x D) matrix, where D = p(p+1)/2
  ### outer loop: iterate over the row of the 2nd order derivative
  ### inner loop: iterate over the col of the 2nd order derivative
  ### We use t_ind in 2 ways:
  ### (1) uniquely identify the row that we are in (1st order derivative)
  ### (2) lets us iterate through the 2nd order derivatives
  ###

  # testblock = matrix(NA, 6, 6)
  stride = D
  block = 1

  for (r in 1:(D-1)) {

    # obtain matrix index of the 1st order derivative (uniquely ID row)
    i = t_ind[r,1] # 1st order row index
    j = t_ind[r,2] # 2nd order col index

    for (c in (r+1):D) { # start on the column after the diagonal

      if (c <= (block * stride)) {
        # obtain matrix index of the 2nd order derivative
        # d^2(psi) / (dpsi_ij dpsi_kl)
        k = t_ind[c,1] # 2nd order row index
        l = t_ind[c,2] # 2nd order col index

        # if (r == 3 && c == 12) {
        #   print(paste('i = ', i, ', j = ', j, ', k = ', k, ', l = ', l, sep = ''))
        # }

        H[r,c] = d2(i, j, k, l, psi_mat, params)
        H[c,r] = H[r,c]
      } else {
        break
      }

    } # end for loop over the columns

    if ((r %% stride) == 0) {
      block = block + 1
    }

    # print(paste("finished row r =", r, sep = ''))
  } # end outer for loop over rows


  return(H)


}






ff = function(u, params) {

  G = params$G
  D = params$D
  t_ind = params$t_ind
  n_nonfree = params$n_nonfree
  vbar = params$vbar
  psi_mat = create_psi_mat(u, params)
  #
  #
  # obj = list(G = G, p = p, psi_mat = psi_mat, n_nonfree = n_nonfree,
  #            vbar = vbar)

  H = matrix(0, D, D)
  for (d in 1:D) {

    i = t_ind[d,1] # row index
    j = t_ind[d,2] # col index

    # populate H[i,i] element
    if (i == j) { # one of the original diag elements in the psi matrix, Psi_ii
      tmp = 0
      for (a in 1:n_nonfree) {
        rr = vbar[a, 1] # row index of non-free element
        ss = vbar[a, 2] # col index of non-free element

        tmp = tmp + d1(rr, ss, i, j, psi_mat, G)^2
      }
      # if (tmp != 0) {print("tmp != 0")}

      H[d,d] = (params$b + params$nu_i[i] - 1) / psi_mat[i,i]^2 + 1 + tmp
    } else { # 2nd order partial for original off-diag elements in the psi mat

      tmp = 0
      for (a in 1:n_nonfree) {
        rr = vbar[a, 1] # row index of non-free element
        ss = vbar[a, 2] # col index of non-free element

        tmp = tmp + d1(rr, ss, i, j, psi_mat, G)^2
      }
      H[d,d] = 1 + tmp
    }
  }


  ### Populate the hessian matrix. This is a (D x D) matrix, where D = p(p+1)/2
  ### outer loop: iterate over the row of the 2nd order derivative
  ### inner loop: iterate over the col of the 2nd order derivative
  ### We use t_ind in 2 ways:
  ### (1) uniquely identify the row that we are in (1st order derivative)
  ### (2) lets us iterate through the 2nd order derivatives
  ###
  for (r in 1:(D-1)) {

    # obtain matrix index of the 1st order derivative (uniquely ID row)
    i = t_ind[r,1] # 1st order row index
    j = t_ind[r,2] # 2nd order col index

    for (c in (r+1):D) { # start on the column after the diagonal
      # obtain matrix index of the 2nd order derivative
      # d^2(psi) / (dpsi_ij dpsi_kl)
      k = t_ind[c,1] # 2nd order row index
      l = t_ind[c,2] # 2nd order col index

      # if (r == 3 && c == 12) {
      #   print(paste('i = ', i, ', j = ', j, ', k = ', k, ', l = ', l, sep = ''))
      # }

      H[r,c] = d2(i, j, k, l, psi_mat, params)
      H[c,r] = H[r,c]
    } # end for loop over the columns
    # print(paste("finished row r =", r, sep = ''))
  } # end outer for loop over rows


  return(H)
}





## this function will call a recursive function
d2 = function(i, j, k, l, psi_mat, params) {

  G = params$G
  tmp = numeric(params$n_nonfree)
  for (n in 1:params$n_nonfree) {
    r = params$vbar[n,1]
    s = params$vbar[n,2]
    if (psi_mat[r,s] == 0) { # avoid needless recursion if coefficient is 0
      tmp[n] = d1(r,s,k,l,psi_mat,G) * d1(r,s,i,j,psi_mat,G)
    } else {
      tmp[n] = d1(r,s,k,l,psi_mat,G) * d1(r,s,i,j,psi_mat,G) +
        psi_mat[r,s] * d2_rs(r,s,i,j,k,l, psi_mat,G)
    }
  } # end for loop iterating over non-free elements
  return(sum(tmp))
} # end d2() function





## function that will recursively compute the 2nd order derivative
## assumption: we don't have i = k AND j = l, i.e., comuting 2nd order deriv.
## for diagonal term in the hessian
d2_rs = function(r, s, i, j, k, l, psi_mat, G) {

  # G = obj$G
  # psi_mat = obj$psi_mat

  if (G[r,s] > 0)                           { return (0) } # free element
  if (r < i || r < k)                       { return (0) } # row below
  if (r == i && r == k && (s < j || s < l)) { return (0) } # same row, col after

  # general case: recursive call for the 2nd order derivative
  tmp = numeric(r-1)
  for (m in 1:(r-1)) {

    if (psi_mat[m,s] == 0 && psi_mat[m,r] == 0) {
      # print('both psi_ms = 0, psi_mr = 0, save recursion')
      tmp[m] = - 1 / psi_mat[r,r] *
        (d1(m,r,i,j,psi_mat,G) * d1(m,s,k,l,psi_mat,G) +
           d1(m,r,k,l,psi_mat,G) * d1(m,s,i,j,psi_mat,G))
    } else {

      if (r == i && i == j) { # case: d^2(psi_rs) / (dpsi_rr dpsi_kl)

        if (psi_mat[m,s] == 0) {
          print("saving ms = 0")
          tmp[m] = 1/psi_mat[r,r]^2 * (psi_mat[m,r] * d1(m,s,k,l,psi_mat,G)) -
            1/psi_mat[r,r] *
            (d2_rs(m,s,i,j,k,l,psi_mat,G) * psi_mat[m,r] +
               d1(m,r,i,j,psi_mat,G) * d1(m,s,k,l,psi_mat,G) +
               d1(m,r,k,l,psi_mat,G) * d1(m,s,i,j,psi_mat,G))
        } else if (psi_mat[m,r] == 0) {
          print("saving mr = 0")
          tmp[m] = 1/psi_mat[r,r]^2 * (d1(m,r,k,l,psi_mat,G) * psi_mat[m,s]) -
            1/psi_mat[r,r] *
            (d2_rs(m,r,i,j,k,l,psi_mat,G) * psi_mat[m,s] +
               d1(m,r,i,j,psi_mat,G) * d1(m,s,k,l,psi_mat,G) +
               d1(m,r,k,l,psi_mat,G) * d1(m,s,i,j,psi_mat,G))
        } else {
          tmp[m] = 1/psi_mat[r,r]^2 *
            (d1(m,r,k,l,psi_mat,G) * psi_mat[m,s] + psi_mat[m,r] * d1(m,s,k,l,psi_mat,G)) -
            1/psi_mat[r,r] *
            (d2_rs(m,r,i,j,k,l,psi_mat,G) * psi_mat[m,s] +
               d2_rs(m,s,i,j,k,l,psi_mat,G) * psi_mat[m,r] +
               d1(m,r,i,j,psi_mat,G) * d1(m,s,k,l,psi_mat,G) +
               d1(m,r,k,l,psi_mat,G) * d1(m,s,i,j,psi_mat,G))
        }

      } else if (r == k && k == l) { # case: d^2(psi_rs) / (dpsi_ij dpsi_rs)
        # flip the order so we get into the first if()
        tmp[m] = d2_rs(r,s,k,l,i,j,psi_mat,G)
      } else {
        ### case: r != i

        if (psi_mat[m,s] == 0) {
          print("saving ms = 0")
          tmp[m] = - 1 / psi_mat[r,r] *
            (d1(m,r,i,j,psi_mat,G) * d1(m,s,k,l,psi_mat,G) +
               d1(m,r,k,l,psi_mat,G) * d1(m,s,i,j,psi_mat,G) +
               psi_mat[m,r] * d2_rs(m,s,i,j,k,l,psi_mat,G))
        } else if (psi_mat[m,r] == 0) {
          print("saving mr = 0")
          tmp[m] = - 1 / psi_mat[r,r] *
            (d1(m,r,i,j,psi_mat,G) * d1(m,s,k,l,psi_mat,G) +
               d1(m,r,k,l,psi_mat,G) * d1(m,s,i,j,psi_mat,G) +
               psi_mat[m,s] * d2_rs(m,r,i,j,k,l,psi_mat,G))
        } else {
          tmp[m] = - 1 / psi_mat[r,r] *
            (d1(m,r,i,j,psi_mat,G) * d1(m,s,k,l,psi_mat,G) +
               d1(m,r,k,l,psi_mat,G) * d1(m,s,i,j,psi_mat,G) +
               psi_mat[m,s] * d2_rs(m,r,i,j,k,l,psi_mat,G) +
               psi_mat[m,r] * d2_rs(m,s,i,j,k,l,psi_mat,G))
        }
      }


    }




  }

  return(sum(tmp))
} # end d2_rs() function










# create_psi_mat(u, params)
# stride = params_G5$D # dimension of parameter in the 1 graph example
# block = 1 # index for which graph we are on
# hess_list = vector("list", length = n_G5)
# inv_hess_list  = vector("list", length = n_G5)
# for (n in 1:n_G5) {
#
#   # extract the first (p1 x p1) block out psi_mat
#   start = (block-1) * stride + 1
#   end   = block * stride
#   u_n = u[start:end]
#   # psi_mat_n = create_psi_mat(u_n, params_G5)
#   # print(psi_mat_n)
#   block = block + 1
#
#   hess_list[[n]]     = ff_fast(u_n, params_G5)
#   inv_hess_list[[n]] = solve(hess_list[[n]])
# }
#
# # return these
# H      = matrix(Matrix::bdiag(hess_list), D, D)
# H_inv  = matrix(Matrix::bdiag(inv_hess_list), D, D)
#
# big_hess      = ff_fast(u, params)
# big_hess_inv  = solve(big_hess)
#
# all.equal(H, big_hess)
# all.equal(H_inv, big_hess_inv)


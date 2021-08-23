

u = u_df[10,1:D_u] %>% unlist %>% unname
#
#
# Rcpp::sourceCpp("C:/Users/ericc/mlike_approx/speedup/hiw.cpp")
#
# cov_loglik(u, params)
# HIW_loglik(u, params)
#
# cov_logprior(u, params)
# HIW_logprior(u, params)
#
# psi(u, params)
#
# u_df_fast = hybridml::preprocess(post_samps, D_u, params)
# u_df_slow = preprocess(post_samps, D_u, params)
#
# microbenchmark(grad(u, params),
#                fast_grad(u, params))
#
# all.equal(grad(u, params),
#                fast_grad(u, params))
#
# all.equal(fast_hess(u, params),
#           hess(u, params))
#
# microbenchmark(hess_cpp = hess(u, params),
#                hess_r = fast_hess(u, params),
#                hess_numer = slow_hess(u, params))
#
#
# Lt = matrix(0, D, D)     # (D x D) lower triangular matrix
# Lt_vec_0 = numeric(D_0)  # (D_0 x 1) vector to fill upper triangular, Lt
# Lt_vec_0[edgeInd] = u
# Lt[upper.tri(Lt, diag = T)] = Lt_vec_0

fast_grad = function(u, params) {

  D   = params$D
  D_0 = params$D_0
  N   = params$N
  S   = params$S

  Lt = matrix(0, D, D)     # (D x D) lower triangular matrix
  Lt_vec_0 = numeric(D_0)  # (D_0 x 1) vector to fill upper triangular, Lt
  Lt_vec_0[edgeInd] = u
  Lt[upper.tri(Lt, diag = T)] = Lt_vec_0

  nu = rowSums(Lt != 0) - 1
  nu[D] = 0
  xi = b + nu - 1

  grad_mat = - diag((xi + N) / diag(Lt)) + Lt + Lt %*% S
  grad_vec = grad_mat[upper.tri(grad_mat, diag = T)]
  # grad_vec[params$edgeInd] # consider only terms that have edge in the graph

  return(grad_vec[params$edgeInd])
}


fast_hess = function(u, params) {

  D   = params$D
  D_0 = params$D_0
  D_u = params$D_u
  S   = params$S

  Lt = matrix(0, D, D)     # (D x D) lower triangular matrix
  Lt_vec_0 = numeric(D_0)  # (D_0 x 1) vector to fill upper triangular, Lt
  Lt_vec_0[edgeInd] = u
  Lt[upper.tri(Lt, diag = T)] = Lt_vec_0

  nu = rowSums(Lt != 0) - 1
  nu[D] = 0
  xi = b + nu - 1

  t_ind = which(Lt != 0, arr.ind = T)

  H = matrix(NA, D_u, D_u)

  for (r in 1:D_u) {

    i = t_ind[r, 1] # row of 1st partial
    j = t_ind[r, 2] # col of 1st partial

    c = r
    while (c <= D_u) {

      k = t_ind[c, 1] # row of 2nd partial
      l = t_ind[c, 2] # col of 2nd partial

      if (i != k) {
        H[r, c] = H[c, r] = 0
      } else if (i == j && k == i && l > j) {
        H[r, c] = H[c, r] = -S[l, j]
      } else if (i == j && j == k && k == l) {
        H[r, c] = H[c, r] = -1/Lt[i,i]^2 * (N + xi[i]) - S[i,i] - 1
      } else if (i != j && k == i && l == j) {
        H[r, c] = H[c, r] = -S[l, j] - 1
      } else if (i != j && k == i && l > j) {
        H[r, c] = H[c, r] = -S[l, j]
      }
      c = c + 1
    }

  }
  return(-H)
}









#
# psi(u, params)
#
# lambda(u, params)
#
#
# Lt = matrix(0, D, D)     # (D x D) lower triangular matrix
# Lt_vec_0 = numeric(D_0)  # (D_0 x 1) vector to fill upper triangular, Lt
# Lt_vec_0[edgeInd] = u
# Lt[upper.tri(Lt, diag = T)] = Lt_vec_0
#
# nu = rowSums(Lt != 0) - 1
# nu[D] = 0
#
# xi = b + nu - 1
#
# grad_mat = - Lt + diag(xi / diag(Lt))
# grad_mat[upper.tri(grad_mat, diag = T)][edgeInd]
#
# pracma::grad(HIW_logprior, u, params = params)
#
#
# # index_mat = matrix(0, D, D)
# # index_mat[upper.tri(index_mat, diag = T)] = 1:D_u
# # index_mat[upper.tri(index_mat, diag = T)]
# # t_ind = which(index_mat!=0,arr.ind = T)
# # t_ind
#
# t_ind = which(Lt != 0, arr.ind = T)
#
# H = matrix(NA, D_u, D_u)
# r = 1
# for (r in 1:D_u) {
#
#     i = t_ind[r, 1] # row of 1st partial
#     j = t_ind[r, 2] # col of 1st partial
#
#     c = r
#     while (c <= D_u) {
#
#         k = t_ind[c, 1] # row of 2nd partial
#         l = t_ind[c, 2] # col of 2nd partial
#
#         if (i != k) {
#           H[r, c] = 0
#         } else if (i == j && k == i && l > j) {
#           H[r, c] = 0
#         } else if (i == j && j == k && k == l) {
#           H[r, c] = - xi[i] / Lt[i, i]^2 - 1
#         } else if (i != j && k == i && l == j) {
#           H[r, c] = -1
#         } else if (i != j && k == i && l > j) {
#           H[r, c] = 0
#         }
#         c = c + 1
#     }
#
# }
#
# H[is.na(H)] = 0
# H1 = Matrix::forceSymmetric(H)
#
#
# pracma::hessian(HIW_logprior, u, params = params) %>% diag
#
# all.equal(pracma::hessian(HIW_logprior, u, params = params), H)
# all.equal(pracma::hessian(HIW_logprior, u, params = params) %>% diag,
#           H %>% diag)
#
# H %>% diag %>% length
#
#
# ### gradient of log likelihood
# # pracma::hessian(HIW_logprior, u, params = params) %>% diag
#
# pracma::grad(HIW_loglik, u, params = params)
#
# dLdLt = diag(N / diag(Lt)) - Lt %*% S
# dLdLt[upper.tri(dLdLt, diag = TRUE)][edgeInd]
#
# all.equal(pracma::grad(HIW_loglik, u, params = params),
#           dLdLt[upper.tri(dLdLt, diag = TRUE)][edgeInd])
#
#
# H_n = pracma::hessian(HIW_loglik, u, params = params)
#
#
# H = matrix(NA, D_u, D_u)
#
# for (r in 1:D_u) {
#
#   i = t_ind[r, 1] # row of 1st partial
#   j = t_ind[r, 2] # col of 1st partial
#
#   c = r
#   while (c <= D_u) {
#
#     k = t_ind[c, 1] # row of 2nd partial
#     l = t_ind[c, 2] # col of 2nd partial
#
#     if (i != k) {
#       H[r, c] = 0
#     } else if (i == j && k == i && l > j) {
#       H[r, c] = -S[l, j]
#     } else if (i == j && j == k && k == l) {
#       H[r, c] = -N/Lt[i,i]^2 - S[i,i]
#     } else if (i != j && k == i && l == j) {
#       H[r, c] = -S[l, j]
#     } else if (i != j && k == i && l > j) {
#       H[r, c] = -S[l, j]
#     }
#     c = c + 1
#   }
#
# }
# H %>% diag
#
# H2 = Matrix::forceSymmetric(H)
#
# all.equal(round(H[upper.tri(H, diag = T)], 4),
#           round(H_n[upper.tri(H_n, diag = T)], 4))
#
# H_0 = - as.matrix(H1 + H2)
# all.equal(H_0,
#           pracma::hessian(psi, u, params = params))
#
#
# # ------------------------------------------------------------------------------
#
# ## gradient of psi:
# pracma::grad(psi, u, params = params)
#
# grad_mat = - diag((xi + N) / diag(Lt)) + Lt + Lt %*% S
# grad_vec = grad_mat[upper.tri(grad_mat, diag = T)]
# grad_vec[edgeInd] # consider only terms that have edge in the graph
#
# all.equal(grad_vec[edgeInd], pracma::grad(psi, u, params = params))
#
#
# ## hessian of psi:
# H = matrix(NA, D_u, D_u)
#
# for (r in 1:D_u) {
#
#   i = t_ind[r, 1] # row of 1st partial
#   j = t_ind[r, 2] # col of 1st partial
#
#   c = r
#   while (c <= D_u) {
#
#     k = t_ind[c, 1] # row of 2nd partial
#     l = t_ind[c, 2] # col of 2nd partial
#
#     if (i != k) {
#       H[r, c] = H[c, r] = 0
#     } else if (i == j && k == i && l > j) {
#       H[r, c] = H[c, r] = -S[l, j]
#     } else if (i == j && j == k && k == l) {
#       H[r, c] = H[c, r] = -1/Lt[i,i]^2 * (N + xi[i]) - S[i,i] - 1
#     } else if (i != j && k == i && l == j) {
#       H[r, c] = H[c, r] = -S[l, j] - 1
#     } else if (i != j && k == i && l > j) {
#       H[r, c] = H[c, r] = -S[l, j]
#     }
#     c = c + 1
#   }
#
# }
# H
# H %>% isSymmetric.matrix()
#
# H_0 = -H
# H_1 = pracma::hessian(psi, u, params = params)
#
# all.equal(-H,
#           H_1)
#
#
#
#
# u = u_df[1,1:D_u] %>% unlist %>% unname
#
# all.equal(grad(u, params),
#           fast_grad(u, params))
#
# library(microbenchmark)
# microbenchmark(grad(u, params),
#                fast_grad(u, params))
#
#
# all.equal(hess(u, params),
#           fast_hess(u, params))
#
# microbenchmark(hess(u, params),
#                fast_hess(u, params),
#                times = 10)

























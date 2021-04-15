


psi = function(u, params) {

  p     = params$p
  G_5   = params$G_5
  b     = params$b
  nu_i  = params$nu_i
  P     = params$P
  b_i   = params$b_i

  ## first reconstruct entire Psi matrix
  Psi_copy = matrix(0, p, p)
  # Phi = matrix(obj$z0, p, p)

  FREE_PARAM_MAT = upper.tri(diag(1, p), diag = T) & G_5
  u_mat = FREE_PARAM_MAT
  u_mat = matrix(0, p, p)
  u_mat[FREE_PARAM_MAT] = u
  u_mat

  ## u_mat should have all the free elements
  for (i in 1:p) {
    for (j in i:p) {
      if (G_5[i,j] > 0) {
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
            tmp2 = u_mat[r,j] + sum(u_mat[r,r:(j-1)] * P[r:(i-1),j] / P[j,j])
            tmp[r] = tmp1 * tmp2
          }
          u_mat[i,j] = x0 - 1 / u_mat[i,i] * sum(tmp)
        }
      }
    }
  }

  t0 = p * log(2) +
    sum((b + b_i - 1) * log(diag(P)) + (b + nu_i - 1) * log(diag(u_mat)))
  t1 = -0.5 * sum(u_mat[upper.tri(u_mat, diag = TRUE)]^2)
  -(t0 + t1)
}


## compute the gwish density
# psi_38 = function(u, params) {
#
#   p     = params$p
#   G_5   = params$G_5
#   b     = params$b
#   nu_i  = params$nu_i
#   P     = params$P
#   b_i   = params$b_i
#
#   ## first reconstruct entire Psi matrix
#   Psi_copy = matrix(0, p, p)
#   # Phi = matrix(obj$z0, p, p)
#
#   FREE_PARAM_MAT = upper.tri(diag(1, p), diag = T) & G_5
#   u_mat = FREE_PARAM_MAT
#   u_mat = matrix(0, p, p)
#   u_mat[FREE_PARAM_MAT] = u
#   u_mat
#
#   ## u_mat should have all the free elements
#   for (i in 1:p) {
#     for (j in i:p) {
#       if (G_5[i,j] > 0) {
#         next # u_mat[i,j] already has value from input
#       } else {
#         if (i == 1) {
#           u_mat[i,j] = -1/P[j,j] * sum(u_mat[i,i:(j-1)] * P[i:(j-1), j])
#         } else {
#           # for rows other than first row
#           x0 = -1/P[j,j] * sum(u_mat[i,i:(j-1)] * P[i:(j-1), j])
#
#           # old formulation that still depends on Phi
#           # Psi_copy[i,j] = x0 -
#           #   sum(Phi[1:(i-1),i] / P[i,i] / (Phi[i,i] / P[j,j]) *
#           #         Phi[1:(i-1), j] / P[j,j])
#
#           # new formulation that uses Eq. (31) in Atay -- no longer uses any
#           # non-free parameters: only uses recursively defined parameters
#           # and elements of the cholesky factorization of the scale matrix
#           tmp = numeric(i-1)
#           for (r in 1:(i-1)) {
#             tmp1 = u_mat[r,i] + sum(u_mat[r,r:(i-1)] * P[r:(i-1),i] / P[i,i])
#             tmp2 = u_mat[r,j] + sum(u_mat[r,r:(j-1)] * P[r:(i-1),j] / P[j,j])
#             tmp[r] = tmp1 * tmp2
#           }
#           u_mat[i,j] = x0 - 1 / u_mat[i,i] * sum(tmp)
#         }
#       }
#     }
#   }
#
#   Psi_copy = u_mat
#
#   # print(u_mat)
#
#   # constant term
#   t0 = sum(0.5 * (b + nu_i) * log(2) + 0.5 * nu_i * log(2 * pi) +
#              lgamma(0.5 * (b + nu_i))) +
#     sum(0.5 * (b + b_i - 1) * log(diag(P)^2))
#
#   ## compute term that sums over non-free terms
#   NON_FREE = !edgeInd
#   UD_FREE  = (UPPER_DIAG & G_5)
#   diag(UD_FREE) = FALSE
#
#   # non-free terms
#   psi_non_free = Psi_copy[UPPER_DIAG][NON_FREE]
#   t1 = - 0.5 * sum(psi_non_free^2)
#
#   # product over diagonal terms
#   t2 = sum(-lgamma(0.5 * (b + nu_i)) +
#              (0.5 * (b + nu_i) - 1) *
#              log(0.5 * diag(Psi_copy)^2) - 0.5 * diag(Psi_copy)^2)
#
#   # product over off-diagonal terms, free terms
#   psi_free_ud = Psi_copy[UD_FREE]
#   t3 = sum(-log(2 * pi) - 0.5 * psi_free_ud^2)
#
#   t_sum = t0 + t1 + t2 + t3
#   return(-t_sum)
#
# }

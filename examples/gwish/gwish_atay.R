
p = 5
G_5 = matrix(c(1,1,0,1,1,
               1,1,1,0,0,
               0,1,1,1,1,
               1,0,1,1,1,
               1,0,1,1,1), p, p)
b = 3
V = diag(1, p)
gnorm(G_5, b, V, 100)

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

u_mat = chol(Omega_G)
u = u_mat[upper.tri(u_mat, diag = TRUE)][edgeInd]

psi(u, params)

pracma::grad(psi, u, params = params)


J = 1000
N = 0
S = matrix(0, p, p)
b = 3
D = sum(edgeInd)

post_gW = sampleGW(J, edgeInd, G_5, b, N, V, S) %>% data.frame()

u_df = preprocess(post_gW, D, params)     # J x (D_u + 1)
u_df %>% head()

u = u_df[1,1:D] %>% unname %>% unlist
u
pracma::hessian(psi, u, params = params)

psi = function(u, params) {

  Phi     = params$Phi
  P       = params$P
  G_5     = params$G_5
  p       = params$p
  edgeInd = params$edgeInd
  b       = params$b
  nu_i    = params$nu_i
  b_i     = params$b_i

  ## form Psi that we can use in the function
  Psi_free = matrix(0, p, p)
  Psi_free[upper.tri(Psi_free, diag = TRUE)][edgeInd] = u

  ## element-wise definition of Psi_{r, s}
  Psi_copy = matrix(0, p, p)
  for (r in 1:p) {
    for (s in r:p) {
      if (G_5[r,s] > 0) {
        Psi_copy[r, s] = Psi_free[r, s]
      }
      if (r == s) {
        Psi_copy[r, s] = Phi[r,s] / P[s, s]
      } else {
        Psi_copy[r, s] = - sum(Psi_copy[r, r:(s-1)] * P[r:(s-1), s] / P[s, s]) +
          Phi[r,s] / P[s, s]
      }
    }
  }

  # constant term
  t0 = sum(0.5 * (b + nu_i) * log(2) + 0.5 * nu_i * log(2 * pi) +
             lgamma(0.5 * (b + nu_i))) +
    sum(0.5 * (b + b_i - 1) * log(diag(P)^2))

  ## compute term that sums over non-free terms
  NON_FREE = !edgeInd
  UD_FREE  = (UPPER_DIAG & G_5)
  diag(UD_FREE) = FALSE

  # non-free terms
  psi_non_free = Psi_copy[UPPER_DIAG][NON_FREE]
  t1 = - 0.5 * sum(psi_non_free^2)

  # product over diagonal terms
  t2 = sum(-lgamma(0.5 * (b + nu_i)) +
             (0.5 * (b + nu_i) - 1) *
             log(0.5 * diag(Psi_copy)^2) - 0.5 * diag(Psi_copy)^2)

  # product over off-diagonal terms, free terms
  psi_free_ud = Psi_copy[UD_FREE]
  t3 = sum(-log(2 * pi) - 0.5 * psi_free_ud^2)

  t_sum = t0 + t1 + t2 + t3

  return(-t_sum)
}






# UPPER_DIAG = upper.tri(Psi_copy, diag = TRUE)
# ## extract the part of Psi that will be used as input into log posterior
# Psi_vec = Psi[upper.tri(Psi, diag = TRUE)]
# Psi_vec = Psi_vec[edgeInd] ## free elements passed into log posterior only
#
# ## form Psi that we can use in the function
# Psi_free = matrix(0, p, p)
# Psi_free[upper.tri(Psi_free, diag = TRUE)][edgeInd] = Psi_vec
# Psi_free
#
# ## element-wise definition of Psi_{r, s}
# Psi_copy = matrix(0, p, p)
# for (r in 1:p) {
#   for (s in r:p) {
#     if (G_5[r,s] > 0) {
#       Psi_copy[r, s] = Psi_free[r, s]
#     }
#     if (r == s) {
#       Psi_copy[r, s] = Phi[r,s] / P[s, s]
#     } else {
#       Psi_copy[r, s] = - sum(Psi_copy[r, r:(s-1)] * P[r:(s-1), s] / P[s, s]) +
#         Phi[r,s] / P[s, s]
#     }
#   }
# }
#
#
# # constant term
# t0 = sum(0.5 * (b + nu_i) * log(2) + 0.5 * nu_i * log(2 * pi) +
#            lgamma(0.5 * (b + nu_i))) +
#   sum(0.5 * (b + b_i - 1) * log(diag(P)^2))
#
# ## compute term that sums over non-free terms
# NON_FREE = !edgeInd
# UD_FREE  = (UPPER_DIAG & G_5)
# diag(UD_FREE) = FALSE
#
#
# # non-free terms
# psi_non_free = Psi_copy[UPPER_DIAG][NON_FREE]
# t1 = - 0.5 * sum(psi_non_free^2)
#
# # product over diagonal terms
# t2 = sum(-lgamma(0.5 * (b + nu_i)) +
#            (0.5 * (b + nu_i) - 1) * log(0.5 * diag(Psi_copy)^2) - 0.5 * diag(Psi_copy)^2)
#
# # product over off-diagonal terms, free terms
# psi_free_ud = Psi_copy[UD_FREE]
# t3 = sum(-log(2 * pi) - 0.5 * psi_free_ud^2)
#
# t0 + t1 + t2 + t3
#
#
# UPPER_DIAG = upper.tri(Psi_copy, diag = TRUE)
#
# ## extract the part of Psi that will be used as input into log posterior
# Psi_vec = Psi[upper.tri(Psi, diag = TRUE)]
# Psi_vec = Psi_vec[edgeInd] ## free elements passed into log posterior only
#

# t1 + t2
# -0.5 * sum(Psi[UPPER_DIAG][NON_FREE]^2) +
#   sum(log(1/gamma((b + nu_i) / 2)) +
#         (0.5 * (b + nu_i) - 1) * log(diag(Psi)^2 / 2) - 0.5 * diag(Psi)^2)

# t3
# sum(log(1/(2*pi)) - 0.5 * Psi[UD_FREE]^2)



all.equal(Psi, Psi_copy)

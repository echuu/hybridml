
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

####

# compute the hessian block by block

u_df %>% head
u = u_df[1,1:D] %>% unlist %>% unname
u

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
  psi_mat_n = create_psi_mat(u_n, params_G5)
  # print(psi_mat_n)
  block = block + 1

  hess_list[[n]]     = ff_fast(u_n, params_G5)
  inv_hess_list[[n]] = solve(hess_list[[n]])
}

H      = matrix(Matrix::bdiag(hess_list), D, D)
H_inv  = matrix(Matrix::bdiag(inv_hess_list), D, D)

big_hess      = ff_fast(u, params)
big_hess_inv  = solve(big_hess)

all.equal(H, big_hess)
all.equal(H_inv, big_hess_inv)


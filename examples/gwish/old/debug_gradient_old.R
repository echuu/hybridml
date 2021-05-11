

r = 2; s = 4;
i = 2; j = 2;
k = 1; l = 5;

1/psi_mat[r,r]^2 * psi_mat[1,r] * psi_mat[1,s] # closed form manual calculation
dpsi_rs(r,s,i,j) # incorrect implementation  = 0
d1(r,s,i,j)      # corrected implementation != 0


dpsi_rs(2,5,2,2)
d1(2,5,2,2)

n_G0 * 12 + n_G1 * 25

## find the indices to extract each block of the hessian that we will invert
## m will iterate over the blocks


u_df %>% head
u = u_df[1,1:D] %>% unlist %>% unname

hess = function(u, params) { pracma::hessian(psi, u, params = params) }
H_numer = hess(u, params)
H_inv_numer = solve(H_numer)

hess = function(u, params) { ff(u, params) }
H = hess(u, params)

H[1:12,13:D]
H[13:D,1:12]

solve(H)[1:12,13:D]
solve(H)[13:D,1:12]

stride = 12
inv_list = vector("list", length = n_G0)
for (m in 1:n_G0) {

  ii = (m-1) * stride + 1
  jj = ii + (stride - 1)
  print(paste("r = ", ii, ", c = ", jj, sep = ''))


  H_m = H_numer[ii:jj, ii:jj]

  inv_list[[m]] = chol2inv(chol(H_m))
}

H_inv = matrix(as.matrix(Matrix::bdiag(inv_list)), D, D)

all.equal(H_inv, H_inv_numer)

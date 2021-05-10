

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

stride = 12
for (m in 1:n_G0) {

  ii = (m-1) * stride + 1
  jj = ii + (stride - 1)
  print(paste("r = ", ii, ", c = ", jj, sep = ''))





}

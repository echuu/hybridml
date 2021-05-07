


i = 3; j = 4;
1 + dpsi_rs(1, 3, i, j)^2 + dpsi_rs(2, 4, i, j)^2 + dpsi_rs(2, 5, i, j)^2



u_df %>% head
u = u_df[1,1:D] %>% unlist %>% unname
psi_mat = create_psi_mat(u, params)
psi_mat

index_mat = matrix(0, p, p)
index_mat[upper.tri(index_mat, diag = T)][edgeInd] = 1:D
index_mat[upper.tri(index_mat, diag = T)]
t_ind = which(index_mat!=0,arr.ind = T)
t_ind

index_mat[lower.tri(index_mat)] = NA
vbar = which(index_mat==0,arr.ind = T) # non-free elements
n_nonfree = nrow(vbar)

if (nrow(t_ind) != D) { print("incorrect # of free parameters")}


## populate the diagonal of the hessian matrix
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

      tmp = tmp + dpsi_rs(rr, ss, i, j)^2
    }
    if (tmp != 0) {print("tmp != 0")}

    H[d,d] = (b + nu_i[i] - 1) / psi_mat[i,i]^2 + 1 + tmp
  } else { # 2nd order partial for original off-diag elements in the psi matrix

    tmp = 0
    for (a in 1:n_nonfree) {
      rr = vbar[a, 1] # row index of non-free element
      ss = vbar[a, 2] # col index of non-free element

      tmp = tmp + dpsi_rs(rr, ss, i, j)^2
    }

    H[d,d] = 1 + tmp
  }
}

data.frame(
  diag(H),
  diag(hess(u, params))
)

round(hess(u, params), 5)


all.equal(diag(H),
          diag(hess(u, params)))



## fill the off-diagonals




hess(u,params)[2,]

i = 1; j = 2;
k = 1; l = 4;
dpsi_rs(1,3,k,l) * dpsi_rs(1,3,i,j) +
  dpsi_rs(2,4,k,l) * dpsi_rs(2,4,i,j) - psi_mat[2,4] / psi_mat[2,2] +
  dpsi_rs(2,5,k,l) * dpsi_rs(2,5,i,j)



i = 1; j = 2;
k = 1; l = 5;
dpsi_rs(1,3,k,l) * dpsi_rs(1,3,i,j) +
  dpsi_rs(2,4,k,l) * dpsi_rs(2,4,i,j) +
  dpsi_rs(2,5,k,l) * dpsi_rs(2,5,i,j) - psi_mat[2,5] / psi_mat[2,2]

d = 1
for (k in 1:(p-1)) {
  for (l in (r+1):(p-1)) {

    i = t_ind[d,1] # row index of first derivative
    j = t_ind[d,2] # col index of first derivative

    # compute d^2(Psi) / d(Psi_ij) d(Psi_kl)
    H[k,l] = d2(i, j, k, l)

  }
}

d2(1,2,1,5)

## this function will call a recursive function
d2 = function(i, j, k, l) {
  tmp = numeric(n_nonfree)
  for (n in 1:n_nonfree) {
    r = vbar[n,1]
    s = vbar[n,2]
    if (psi_mat[r,s] == 0) { # avoid needless recursion if coefficient is 0
      tmp[n] = dpsi_rs(r,s,k,l) * dpsi_rs(r,s,i,j)
    } else {
      tmp[n] = dpsi_rs(r,s,k,l) * dpsi_rs(r,s,i,j) + psi_mat[r,s] * d2_rs(r,s,i,j,k,l)
    }
  } # end for loop iterating over non-free elements
  return(sum(tmp))
} # end d2() function


## function that will recursively compute the 2nd order derivative
## assumption: we don't have i = k AND j = l, i.e., comuting 2nd order deriv.
## for diagonal term in the hessian
d2_rs = function(r, s, i, j, k, l) {

  if (G[r,s] > 0)                           { return (0) } # free element
  if (r < i || r < k)                       { return (0) } # row below
  if (r == i && r == k && (s < j || s < l)) { return (0) } # same row, col after

  # general case: recursive call for the 2nd order derivative
  tmp = numeric(r-1)
  for (m in 1:(r-1)) {
    tmp[m] = - 1 / psi_mat[r,r] *
      dpsi_rs(m,r,i,j) * dpsi_rs(m,s,k,l) + dpsi_rs(m,r,k,l) * dpsi_rs(m,s,i,j) +
      d2_rs(m,r,i,j,k,l) + d2_rs(m,s,i,j,k,l)
  }

  return(sum(tmp))
} # end d2_rs() function

# d^2(psi_25) / (dpsi_12 dpsi_15) = - 1/psi_22
-1/psi_mat[2,2]
r = 2; s = 4; i = 1; j = 2; k = 1; l = 4;
d2_rs(r,s,i,j,k,l)




if (G[r, s] > 0) {
  if (r == i && s == j) { #
    return(0)
  } else {                #
    return(0)
  }
}

if (i > r)                            { return(0) } # (i,j) comes after (r,s)
if (i == r && j > s)                  { return(0) } # same row, but after
if (i == r && j == s && G[r, s] == 0) { return(0) } # deriv wrt to non-free: 0
if (r == 1 && s > r)                  { return(0) } # first row case



if (r > 1) {




}

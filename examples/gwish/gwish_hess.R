


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

      tmp = tmp + d1(obj, rr, ss, i, j)^2
    }
    if (tmp != 0) {print("tmp != 0")}

    H[d,d] = (b + nu_i[i] - 1) / psi_mat[i,i]^2 + 1 + tmp
  } else { # 2nd order partial for original off-diag elements in the psi matrix

    tmp = 0
    for (a in 1:n_nonfree) {
      rr = vbar[a, 1] # row index of non-free element
      ss = vbar[a, 2] # col index of non-free element

      tmp = tmp + d1(obj, rr, ss, i, j)^2
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

d2(i,j,k,l)


i = 1; j = 2;
k = 1; l = 5;
dpsi_rs(1,3,k,l) * dpsi_rs(1,3,i,j) +
  dpsi_rs(2,4,k,l) * dpsi_rs(2,4,i,j) +
  dpsi_rs(2,5,k,l) * dpsi_rs(2,5,i,j) - psi_mat[2,5] / psi_mat[2,2]


#### implementation starts here ------------------------------------------------

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

    H[r,c] = d2(i, j, k, l)
    H[c,r] = H[r,c]
  } # end for loop over the columns
} # end outer for loop over rows

H
round(hess(u, params), 6)

d2(1,2,1,5)

## this function will call a recursive function
d2 = function(i, j, k, l) {
  tmp = numeric(n_nonfree)
  for (n in 1:n_nonfree) {
    r = vbar[n,1]
    s = vbar[n,2]
    if (psi_mat[r,s] == 0) { # avoid needless recursion if coefficient is 0
      tmp[n] = d1(r,s,k,l) * d1(r,s,i,j)
    } else {
      tmp[n] = d1(r,s,k,l) * d1(r,s,i,j) + psi_mat[r,s] * d2_rs(r,s,i,j,k,l)
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

    if (r == i && i == j) { # case: d^2(psi_rs) / (dpsi_rr dpsi_kl)
      tmp[m] = 1/psi_mat[r,r]^2 *
        (d1(m,r,k,l) * psi_mat[m,s] + psi_mat[m,r] * d1(m,s,k,l)) -
        1/psi_mat[r,r] *
        (d2_rs(m,r,i,j,k,l) * psi_mat[m,s] + d2_rs(m,s,i,j,k,l) * psi_mat[m,r] +
            d1(m,r,i,j) * d1(m,s,k,l) + d1(m,r,k,l) * d1(m,s,i,j))
    } else if (r == k && k == l) { # case: d^2(psi_rs) / (dpsi_ij dpsi_rs)
      tmp[m] = d2_rs(r,s,k,l,i,j) # flip the order so we get into the first if()
    } else {
      ### case: r != i
      tmp[m] = - 1 / psi_mat[r,r] *
        (d1(m,r,i,j) * d1(m,s,k,l) + d1(m,r,k,l) * d1(m,s,i,j) +
           d2_rs(m,r,i,j,k,l) + d2_rs(m,s,i,j,k,l))
    }

  }

  return(sum(tmp))
} # end d2_rs() function



round(hess(u, params), 6)
## manual calculation for: d psi / ( d(psi_12) d(psi_22) )
- psi_mat[1,2] * psi_mat[1,4]^2 / psi_mat[2,2]^3 +
  psi_mat[2,4] * psi_mat[1,4] / psi_mat[2,2]^2 -
  psi_mat[1,2] * psi_mat[1,5]^2 / psi_mat[2,2]^3 +
  psi_mat[2,5] * psi_mat[1,5] / psi_mat[2,2]^2
# -0.00006588557

# compare to closed form calculation:
i = 1; j = 2; k = 2; l = 2;
d2(i, j, k, l)

## check closed form
r = 2
1/psi_mat[r,r]^2 * psi_mat[1,2] * psi_mat[2,4] +
  1/psi_mat[2,2]^2 * psi_mat[1,2] * psi_mat[1,4]


H[3,5]
hess(u, params)[3,5]


## another check dpsi / dpsi_22 dpsi_14
H[3,6]
hess(u, params)[3,6]
-psi_mat[1,2]^2 * psi_mat[1,4] / psi_mat[2,2]^3 + psi_mat[2,4] * psi_mat[1,2] / psi_mat[2,2]^2


## another check dpsi / dpsi_22 dpsi_15
H[3,9]
hess(u, params)[3,9]
-psi_mat[1,2]^2 * psi_mat[1,5] / psi_mat[2,2]^3 + psi_mat[2,5] * psi_mat[1,2] / psi_mat[2,2]^2

## another check dpsi / dpsi_22 dpsi_55
H[3,12]
hess(u, params)[3,12]

H_0 = hess(u, params)
all.equal(
  H_0[upper.tri(H_0)],
  H[upper.tri(H)]
)


## matching terms:
# (1)
d1(2,4,1,2) * d1(2,4,2,2)
-psi_mat[1,4]^2 * psi_mat[1,2] / psi_mat[2,2]^3

# (3)
d1(2,5,1,2) * d1(2,5,2,2)
-psi_mat[1,5]^2 * psi_mat[1,2] / psi_mat[2,2]^3

# (2) currently returning the wrong thing
psi_mat[2,4] * d2_rs(2,4,2,2,1,2)

r = 2; s = 4; i = 2; j = 2; k = 1; l = 2;
r = 2; s = 4; i = 1; j = 2; k = 2; l = 2;
r == i
# closed form of d^2(psi_24) / (dpsi_22 dpsi_12)
d2_rs(r,s,i,j,k,l)

## true value of d^2(psi_24) / (dpsi_22 dpsi_12)
psi_mat[1,4] / psi_mat[2,2]^2

## check this one because it will get flagged by the first if, but it technically
## belongs in the last else()
r = 2; s = 4; i = 2; j = 2; k = 1; l = 3;

















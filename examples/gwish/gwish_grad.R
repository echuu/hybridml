

# gradient and hessian for gwish()

Psi_mat = matrix(0, p, p)

Psi_mat[upper.tri(Psi_mat, diag = TRUE)][edgeInd] = u_k

Psi_mat & G_5


u_df %>% head
u = u_df[1,1:D] %>% unlist %>% unname
psi_mat = create_psi_mat(u, params)
psi_mat

i = 1; j = 1;
(b + nu_i[i] - 1) / psi_mat[i,i]^2 + 1


diag(hess(u, params))

create_psi_mat(u, params)
matrix(samps$Psi[1,], p, p)

grad(u, params)


P = chol(solve(V)) # upper cholesky factor; D^(-1) = TT'  in Atay paper

# used to find index of elements
index_mat = matrix(0, p, p)
index_mat[upper.tri(index_mat, diag = T)][edgeInd] = 1:D
index_mat[upper.tri(index_mat, diag = T)]
t_ind = which(index_mat!=0,arr.ind = T)
t_ind
# these are the (row, col) index for the free-parameters
# use these to determine which 2nd order partials to diff wrt


# h = function(i, j) {
#   return(P[i,j] / P[j, j])
# }
#
# psi_mat = create_psi_mat(u, params)
#
# psi_mat[1,3] * (- P[1, 3] / P[3, 3])
# - psi_mat[2,4] / psi_mat[2, 2] *
#   (P[1,2] / P[2,2] * (psi_mat[1,4] + sum(psi_mat[1,1:3] * P[1:3,4] / P[4,4])) +
#      (psi_mat[1, 2] + psi_mat[1,1] * P[1,2]/P[2,2]) * P[3, 4] / P[4,4] * (- P[1,3] / P[3,3]))
# psi_mat[2,5]
#
#
#
# (b + nu_i[1] - 1) / psi_mat[1, 1] - psi_mat[1,1]
# (b + nu_i[2] - 1) / psi_mat[2, 2] - psi_mat[2,2]
# (b + nu_i[3] - 1) / psi_mat[3, 3] - psi_mat[3,3]
#
# psi_mat[2,4] * (-1/psi_mat[2,2] * psi_mat[1,4]) - 1/psi_mat[2,2] * psi_mat[1,5] * psi_mat[2,5]
#
# psi_mat[1,2] - psi_mat[2,4] * psi_mat[1,4] / psi_mat[2,2] -
#   psi_mat[2,5] * psi_mat[1,5] / psi_mat[2,2]


## convert gradient matrix into matrix form for easier comparison
tmp = matrix(NA, p, p)
tmp[upper.tri(tmp, diag = TRUE)][edgeInd] = grad(u, params)
tmp

hess(u, params)


psi_mat = create_psi_mat(u, params)
all.equal(psi_mat, matrix(samps$Psi[1,], p, p))

# G = G_5

f = function() {
  gg = matrix(0, p, p)
  for (i in 1:p) {
    for (j in i:p) {
      if (G[i, j] == 0) {
        gg[i,j] = 0
      } else {
        gg[i,j] = dpsi(i,j) ## dpsi / dpsi_ij
      }
    }
  }

  return(gg[upper.tri(gg, diag = TRUE)][edgeInd])
}


library(microbenchmark)
## closed form is about 100x faster, closed form can def be made faster too
microbenchmark(f1 = f(), f2 = grad(u, params))
gg
tmp = matrix(0, p, p)
tmp[upper.tri(tmp, diag = TRUE)][edgeInd] = grad(u, params)
tmp

## extract only free elements into gradient vector
gg[upper.tri(gg, diag = TRUE)][edgeInd]
tmp[upper.tri(tmp, diag = TRUE)][edgeInd]

data.frame(
  form = gg[upper.tri(gg, diag = TRUE)][edgeInd],
  numer = tmp[upper.tri(tmp, diag = TRUE)][edgeInd]
) %>% mutate(diff = abs(form - numer))



all.equal(gg, tmp)

dpsi(1, 1)
i = 2; j = 4;
i = 1; j = 5;
dpsi(i, j)


grad(u, params)

i = 1; j = 5;
dpsi(i, j)


dpsi(1,1)


dpsi = function(i, j) {
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
          d_ij = d_ij + psi_mat[r,s] * dpsi_rs(r, s, i, j)
        }
      }
    }

    return(d_ij - (b + nu_i[i] - 1) / psi_mat[i, i] + psi_mat[i, i])
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
          d_ij = d_ij + psi_mat[r,s] * dpsi_rs(r, s, i, j)
        }
      } # end loop over s
    } # end loop over r

  } # end else

  return((d_ij + psi_mat[i,j]))
}


## compute the derivative: d psi_{rs} / d psi_{ij}
dpsi_rs = function(r, s, i, j) {

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

      tmp_sum[k] = dpsi_rs(k, r, i, j) * psi_mat[k, s] +
        dpsi_rs(k, s, i, j) * psi_mat[k, r]
    }

    return(-1/psi_mat[r,r] * sum(tmp_sum)) ## expression derived from Eq. (99)

  } else {
    print("case not accounted for")
    return(NA)
  }

}




grad_vec = function(u, params) {




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





























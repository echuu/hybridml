

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






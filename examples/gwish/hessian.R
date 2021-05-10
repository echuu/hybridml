

# index_mat = matrix(0, p, p)
# index_mat[upper.tri(index_mat, diag = T)][edgeInd] = 1:D
# index_mat[upper.tri(index_mat, diag = T)]
# t_ind = which(index_mat!=0,arr.ind = T)
# t_ind
#
# index_mat[lower.tri(index_mat)] = NA
# vbar = which(index_mat==0,arr.ind = T) # non-free elements
# n_nonfree = nrow(vbar)


# params = list(G = G_5, P = P, p = p, D = D, edgeInd = edgeInd,
#               b = b, nu_i = nu_i, b_i = b_i,
#               n_nonfree = n_nonfree, vbar = vbar, t_ind = t_ind)
#
#
# hess(u, params)
# ff(u, params)

# install.packages("profvis")
# library(profvis)
# profvis(ff(u, params))



testblock = matrix(0, 6, 6)
D = ncol(testblock)
stride = 2
block = 1
for (r in 1:(D-1)) {

  for (c in (r+1):D) {
    if (c <= (block * stride)) {
      testblock[r,c] = NA
      testblock[c,r] = NA
    } else {
      break
    }
  }

  if ((r %% stride) == 0) {
    block = block + 1
  }

}
testblock



ff_fast = function(u, params) {

  G = params$G
  D = params$D
  t_ind = params$t_ind
  n_nonfree = params$n_nonfree
  vbar = params$vbar
  psi_mat = create_psi_mat(u, params)
  #
  #
  # obj = list(G = G, p = p, psi_mat = psi_mat, n_nonfree = n_nonfree,
  #            vbar = vbar)

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

        tmp = tmp + d1(rr, ss, i, j, psi_mat, G)^2
      }
      # if (tmp != 0) {print("tmp != 0")}

      H[d,d] = (params$b + params$nu_i[i] - 1) / psi_mat[i,i]^2 + 1 + tmp
    } else { # 2nd order partial for original off-diag elements in the psi mat

      tmp = 0
      for (a in 1:n_nonfree) {
        rr = vbar[a, 1] # row index of non-free element
        ss = vbar[a, 2] # col index of non-free element

        tmp = tmp + d1(rr, ss, i, j, psi_mat, G)^2
      }
      H[d,d] = 1 + tmp
    }
  }


  ### Populate the hessian matrix. This is a (D x D) matrix, where D = p(p+1)/2
  ### outer loop: iterate over the row of the 2nd order derivative
  ### inner loop: iterate over the col of the 2nd order derivative
  ### We use t_ind in 2 ways:
  ### (1) uniquely identify the row that we are in (1st order derivative)
  ### (2) lets us iterate through the 2nd order derivatives
  ###

  # testblock = matrix(NA, 6, 6)
  stride = D
  block = 1

  for (r in 1:(D-1)) {

    # obtain matrix index of the 1st order derivative (uniquely ID row)
    i = t_ind[r,1] # 1st order row index
    j = t_ind[r,2] # 2nd order col index

    for (c in (r+1):D) { # start on the column after the diagonal

      if (c <= (block * stride)) {
        # obtain matrix index of the 2nd order derivative
        # d^2(psi) / (dpsi_ij dpsi_kl)
        k = t_ind[c,1] # 2nd order row index
        l = t_ind[c,2] # 2nd order col index

        # if (r == 3 && c == 12) {
        #   print(paste('i = ', i, ', j = ', j, ', k = ', k, ', l = ', l, sep = ''))
        # }

        H[r,c] = d2(i, j, k, l, psi_mat, params)
        H[c,r] = H[r,c]
      } else {
        break
      }

    } # end for loop over the columns

    if ((r %% stride) == 0) {
      block = block + 1
    }

    # print(paste("finished row r =", r, sep = ''))
  } # end outer for loop over rows


  return(H)


}






ff = function(u, params) {

  G = params$G
  D = params$D
  t_ind = params$t_ind
  n_nonfree = params$n_nonfree
  vbar = params$vbar
  psi_mat = create_psi_mat(u, params)
  #
  #
  # obj = list(G = G, p = p, psi_mat = psi_mat, n_nonfree = n_nonfree,
  #            vbar = vbar)

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

        tmp = tmp + d1(rr, ss, i, j, psi_mat, G)^2
      }
      # if (tmp != 0) {print("tmp != 0")}

      H[d,d] = (params$b + params$nu_i[i] - 1) / psi_mat[i,i]^2 + 1 + tmp
    } else { # 2nd order partial for original off-diag elements in the psi mat

      tmp = 0
      for (a in 1:n_nonfree) {
        rr = vbar[a, 1] # row index of non-free element
        ss = vbar[a, 2] # col index of non-free element

        tmp = tmp + d1(rr, ss, i, j, psi_mat, G)^2
      }
      H[d,d] = 1 + tmp
    }
  }


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

      H[r,c] = d2(i, j, k, l, psi_mat, params)
      H[c,r] = H[r,c]
    } # end for loop over the columns
    # print(paste("finished row r =", r, sep = ''))
  } # end outer for loop over rows


  return(H)
}





## this function will call a recursive function
d2 = function(i, j, k, l, psi_mat, params) {

  G = params$G
  tmp = numeric(params$n_nonfree)
  for (n in 1:params$n_nonfree) {
    r = params$vbar[n,1]
    s = params$vbar[n,2]
    if (psi_mat[r,s] == 0) { # avoid needless recursion if coefficient is 0
      tmp[n] = d1(r,s,k,l,psi_mat,G) * d1(r,s,i,j,psi_mat,G)
    } else {
      tmp[n] = d1(r,s,k,l,psi_mat,G) * d1(r,s,i,j,psi_mat,G) +
        psi_mat[r,s] * d2_rs(r,s,i,j,k,l, psi_mat,G)
    }
  } # end for loop iterating over non-free elements
  return(sum(tmp))
} # end d2() function





## function that will recursively compute the 2nd order derivative
## assumption: we don't have i = k AND j = l, i.e., comuting 2nd order deriv.
## for diagonal term in the hessian
d2_rs = function(r, s, i, j, k, l, psi_mat, G) {

  # G = obj$G
  # psi_mat = obj$psi_mat

  if (G[r,s] > 0)                           { return (0) } # free element
  if (r < i || r < k)                       { return (0) } # row below
  if (r == i && r == k && (s < j || s < l)) { return (0) } # same row, col after

  # general case: recursive call for the 2nd order derivative
  tmp = numeric(r-1)
  for (m in 1:(r-1)) {

    if (psi_mat[m,s] == 0 && psi_mat[m,r] == 0) {
      # print('both psi_ms = 0, psi_mr = 0, save recursion')
      tmp[m] = - 1 / psi_mat[r,r] *
        (d1(m,r,i,j,psi_mat,G) * d1(m,s,k,l,psi_mat,G) +
           d1(m,r,k,l,psi_mat,G) * d1(m,s,i,j,psi_mat,G))
    } else {

      if (r == i && i == j) { # case: d^2(psi_rs) / (dpsi_rr dpsi_kl)

        if (psi_mat[m,s] == 0) {
          print("saving ms = 0")
          tmp[m] = 1/psi_mat[r,r]^2 * (psi_mat[m,r] * d1(m,s,k,l,psi_mat,G)) -
            1/psi_mat[r,r] *
            (d2_rs(m,s,i,j,k,l,psi_mat,G) * psi_mat[m,r] +
               d1(m,r,i,j,psi_mat,G) * d1(m,s,k,l,psi_mat,G) +
               d1(m,r,k,l,psi_mat,G) * d1(m,s,i,j,psi_mat,G))
        } else if (psi_mat[m,r] == 0) {
          print("saving mr = 0")
          tmp[m] = 1/psi_mat[r,r]^2 * (d1(m,r,k,l,psi_mat,G) * psi_mat[m,s]) -
            1/psi_mat[r,r] *
            (d2_rs(m,r,i,j,k,l,psi_mat,G) * psi_mat[m,s] +
               d1(m,r,i,j,psi_mat,G) * d1(m,s,k,l,psi_mat,G) +
               d1(m,r,k,l,psi_mat,G) * d1(m,s,i,j,psi_mat,G))
        } else {
          tmp[m] = 1/psi_mat[r,r]^2 *
            (d1(m,r,k,l,psi_mat,G) * psi_mat[m,s] + psi_mat[m,r] * d1(m,s,k,l,psi_mat,G)) -
            1/psi_mat[r,r] *
            (d2_rs(m,r,i,j,k,l,psi_mat,G) * psi_mat[m,s] +
               d2_rs(m,s,i,j,k,l,psi_mat,G) * psi_mat[m,r] +
               d1(m,r,i,j,psi_mat,G) * d1(m,s,k,l,psi_mat,G) +
               d1(m,r,k,l,psi_mat,G) * d1(m,s,i,j,psi_mat,G))
        }

      } else if (r == k && k == l) { # case: d^2(psi_rs) / (dpsi_ij dpsi_rs)
        tmp[m] = d2_rs(r,s,k,l,i,j,psi_mat,G) # flip the order so we get into the first if()
      } else {
        ### case: r != i

        if (psi_mat[m,s] == 0) {
          print("saving ms = 0")
          tmp[m] = - 1 / psi_mat[r,r] *
            (d1(m,r,i,j,psi_mat,G) * d1(m,s,k,l,psi_mat,G) +
               d1(m,r,k,l,psi_mat,G) * d1(m,s,i,j,psi_mat,G) +
               psi_mat[m,r] * d2_rs(m,s,i,j,k,l,psi_mat,G))
        } else if (psi_mat[m,r] == 0) {
          print("saving mr = 0")
          tmp[m] = - 1 / psi_mat[r,r] *
            (d1(m,r,i,j,psi_mat,G) * d1(m,s,k,l,psi_mat,G) +
               d1(m,r,k,l,psi_mat,G) * d1(m,s,i,j,psi_mat,G) +
               psi_mat[m,s] * d2_rs(m,r,i,j,k,l,psi_mat,G))
        } else {
          tmp[m] = - 1 / psi_mat[r,r] *
            (d1(m,r,i,j,psi_mat,G) * d1(m,s,k,l,psi_mat,G) +
               d1(m,r,k,l,psi_mat,G) * d1(m,s,i,j,psi_mat,G) +
               psi_mat[m,s] * d2_rs(m,r,i,j,k,l,psi_mat,G) +
               psi_mat[m,r] * d2_rs(m,s,i,j,k,l,psi_mat,G))
        }
      }


    }




  }

  return(sum(tmp))
} # end d2_rs() function





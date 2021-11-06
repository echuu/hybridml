


psi = function(u, params) {

  p     = params$p
  G     = params$G
  b     = params$b
  nu_i  = params$nu_i
  P     = params$P
  b_i   = params$b_i

  # print(paste("p = ", p, sep = ''))
  # print(paste("b = ", b, sep = ''))
  # cat("nu_i =", nu_i, '\n')
  # cat("P = ", '\n')
  # print(P)
  # cat("b_i =", b_i, '\n')


  ## first reconstruct entire Psi matrix
  # Psi_copy = matrix(0, p, p)
  # Phi = matrix(obj$z0, p, p)

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

  # print(u_mat)

  t0 = p * log(2) +
    sum((b + b_i - 1) * log(diag(P)) + (b + nu_i - 1) * log(diag(u_mat)))
  t1 = -0.5 * sum(u_mat[upper.tri(u_mat, diag = TRUE)]^2)
  -(t0 + t1)
}




# log multivariate gamma function Gamma_p(a)
logmultigamma = function(p, a){
  f = 0.25*p*(p-1)*log(pi)
  for(i in 1:p){ f = f + lgamma(a+0.5-0.5*i) }
  return(f)
}



# extract only the free elements from Psi in vector form
extractFree = function(u, free_ind) {
  u[free_ind]
}



# compute Psi = Phi * P^(-1)
computePsi = function(phi_vec, P_inv) {
  p = nrow(P_inv)
  Phi = matrix(phi_vec, p, p, byrow = FALSE)
  Phi %*% P_inv
}




## sampleGW function: returns the following:
## (1) Psi: free parameters stored as a row vector; this is fed into hybrid algo
## (2) Phi: entire (p x p) matrix where Phi'Phi = Omega, Omega ~ GW(b, V)
samplegw = function(J, G, b, N, V, S, P_inv, param_ind) {

  ## sample from gw distribution: stored as J - (D x D) matrices
  Omega_post = rgwish(J, G, b + N, V + S)

  ## Compute Phi (upper triangular), stored column-wise, Phi'Phi = Omega
  Phi = apply(Omega_post, 3, chol) # (p^2) x J

  ## Compute Psi
  Psi = apply(Phi, 2, computePsi, P_inv = P_inv)

  ## Extract free elements of Psi
  Psi_free = apply(Psi, 2, extractFree, free_ind = param_ind)

  out = list(Phi = t(Phi),
             Psi = t(Psi),
             Psi_free = t(Psi_free),
             Omega = Omega_post)

}


# sampleParams = function(J, G, b, N, V, S) {
#     # array, where each element is a (p x p) precision matrix
#     Omega_post = rgwish(J, G, b + N, V + S) # J x (D x D)
#     post_samps = t(apply(Omega_post, 3, process_samps, edgeIndex = edgeIndex))
# }
#
#
# process_samps = function(Omega, edgeIndex){
#     Lt = chol(Omega)
#     Lt_vec = Lt[upper.tri(Lt, diag = T)]
#     Lt_vec[edgeIndex]
# }
#
#
# sampleGW = function(J, edgeIndex, G, b, N, V, S) {
#
#     Omega_post    = vector("list", J) # store posterior samples in matrix form
#     # Lt_post       = vector("list", J) # store lower cholesky factor
#     # post_samps_0  = matrix(0, J, D_0) # store ENTIRE upper diag in vector form
#     # post_samps    = matrix(0, J, D_u) # store NONZERO upper diag in vector form
#
#     Omega_post = rgwish(J, G, b + N, V + S) # J x (D x D)
#     post_samps = t(apply(Omega_post, 3, process_samps, edgeIndex = edgeIndex))
#     return(post_samps)
# }



preprocess = function(post_samps, D, params) {

  psi_u = apply(post_samps, 1, psi, params = params) %>% unname() # (J x 1)

  # (1.2) name columns so that values can be extracted by partition.R
  u_df_names = character(D + 1)
  for (d in 1:D) {
    u_df_names[d] = paste("u", d, sep = '')
  }
  u_df_names[D + 1] = "psi_u"

  # populate u_df
  u_df = cbind(post_samps, psi_u) # J x (D + 1)
  names(u_df) = u_df_names


  return(u_df)

} # end of preprocess() function -----------------------------------------------




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



# log multivariate gamma function Gamma_p(a)
logmultigamma = function(p, a){
  f = 0.25*p*(p-1)*log(pi)
  for(i in 1:p){ f = f + lgamma(a+0.5-0.5*i) }
  return(f)
}

logfrac = function(dim, b, D){
  temp = b+dim-1
  logfrac = 0.5*temp*log(det(D)) - 0.5*temp*dim*log(2) -
    logmultigamma(dim, 0.5*temp)
  return(logfrac)
}

# log normalizing constant for HIW
logHIWnorm = function(G, b, D){
  junct = makedecompgraph(G)
  cliques = junct$C; separators = junct$S
  nc = length(cliques); ns = length(separators)

  Cnorm = 0
  for(i in 1:nc){
    ind = cliques[[i]]
    Cnorm = Cnorm + logfrac(length(ind), b, D[ind, ind, drop = FALSE])
  }

  Snorm = 0
  if(ns>1){
    for(i in 2:ns){
      ind = separators[[i]]
      Snorm = Snorm + logfrac(length(ind), b, D[ind, ind, drop = FALSE])
    }
  }

  logHIWnorm = Cnorm - Snorm
  return(logHIWnorm)
}

# # log marginal likelihood log(f(Y|G))
# logmarginal = function(Y, G, b, D){
#     n = nrow(Y); p = ncol(Y); S = t(Y)%*%Y
#     logmarginal = -0.5*n*p*log(2*pi) + logHIWnorm(G, b, D) -
#         logHIWnorm(G, b+n, D+S)
#     return(logmarginal)
# }

logmarginal = function(Y, G, b, D, S){
  n = nrow(Y); p = ncol(Y); # S = t(Y)%*%Y
  logmarginal = -0.5*n*p*log(2*pi) + logHIWnorm(G, b, D) -
    logHIWnorm(G, b+n, D+S)
  return(logmarginal)
}



# inverse Wishart density
logSingle = function(dim, b, D, Sigma){
  temp = b + dim - 1
  logSingle = 0.5 * temp * log(det(D)) - 0.5 * temp * dim * log(2) -
    logmultigamma(dim, 0.5 * temp)
  logSingle = logSingle - 0.5 * (b + 2 * dim) * log(det(Sigma)) -
    0.5 * sum(diag( solve(Sigma) %*% D ))

  return(logSingle)
}

# HIW density
logHIW = function(G, b, D, Sigma){
  junct = makedecompgraph(G)
  cliques = junct$C; separators = junct$S
  nc = length(cliques); ns = length(separators)

  Cnorm = 0
  for(i in 1:nc){
    ind = cliques[[i]]
    Cnorm = Cnorm + logSingle(length(ind), b,
                              D[ind, ind, drop = FALSE],
                              Sigma[ind, ind, drop = FALSE])
  }

  Snorm = 0
  if(ns>1){
    for(i in 2:ns){
      ind = separators[[i]]
      Snorm = Snorm + logSingle(length(ind), b,
                                D[ind, ind, drop = FALSE],
                                Sigma[ind, ind, drop = FALSE])
    }
  }

  logHIW = Cnorm - Snorm
  return(logHIW)
}


# Sample the HIW_G(bG,DG) distribution on a graph G with adjacency matrix Adj
HIWsim = function(Adj, bG, DG){
  # check if "Adj" is a matrix object
  if(is.matrix(Adj)==FALSE) { stop("Adj must be a matrix object!") }
  # check if "Adj" is a square matrix
  if(dim(Adj)[1]!=dim(Adj)[2]) { stop("Adj must be a square matrix") }
  # check if "Adj" is symmetric
  if(isSymmetric.matrix(Adj)==FALSE) { stop("Adj must be a symmetric matrix") }

  # check if "DG" is a matrix object
  if(is.matrix(DG)==FALSE) { stop("DG must be a matrix object!") }
  # check if "DG" is a square matrix
  if(dim(DG)[1]!=dim(DG)[2]) { stop("DG must be a square matrix") }
  # check if "DG" is symmetric
  if(isSymmetric.matrix(DG)==FALSE) { stop("DG must be a symmetric matrix") }

  # check if "bG" is greater than 2
  if(bG<=2) { stop("bG must be greater than 2") }
  # check if "Adj" and "DG" are the same size
  if(nrow(Adj)!=nrow(DG)) { stop("Adj and DG must have the same dimension") }
  rMNorm = function(m, V){ return(m+t(chol(V))%*%rnorm(length(m))) }

  p = nrow(Adj)
  temp = makedecompgraph(Adj)
  Cliques = temp$C
  Separators = temp$S
  numberofcliques = length(Cliques)

  ############################################################
  # Creat some working arrays that are computed only once
  C1 = solve(DG[Cliques[[1]], Cliques[[1]]]/bG)
  c1 = Cliques[[1]]
  UN = c1
  DSi = DRS = mU = list()

  for(i in 2:numberofcliques){
    sid = Separators[[i]]
    DSi[[i]] = solve(DG[sid, sid])
    cid = Cliques[[i]]
    dif = sort(setdiff(cid, UN))
    UN = sort(union(cid, UN)) # no need to sort, just playing safe
    sizedif = length(dif)
    DRS[[i]] = DG[dif, dif] - DG[dif, sid] %*% DSi[[i]] %*% DG[sid, dif]
    DRS[[i]] = ( DRS[[i]] + t(DRS[[i]]) )/2
    mU[[i]] = DG[dif, sid] %*% DSi[[i]]
  }

  ############################################################
  # MC Sampling
  UN = c1
  Sigmaj = matrix(0, p, p)
  # sample variance mx on first component
  Sigmaj[c1, c1] = solve(Wishart_InvA_RNG( bG+length(Cliques[[1]])-1, DG[Cliques[[1]], Cliques[[1]]] ))

  for(i in 2:numberofcliques){ # visit components and separators in turn
    dif = sort(setdiff(Cliques[[i]], UN))
    UN = sort(union(Cliques[[i]], UN)) # probably no need to sort, just playing safe
    sizedif = length(dif)
    sid = Separators[[i]]
    SigRS = solve(Wishart_InvA_RNG( bG+length(Cliques[[i]])-1, DRS[[i]] ))
    Ui = rMNorm( as.vector(t(mU[[i]])), kronecker(SigRS, DSi[[i]]))
    Sigmaj[dif, sid] = t(matrix(Ui, ncol = sizedif)) %*% Sigmaj[sid, sid]
    Sigmaj[sid, dif] = t(Sigmaj[dif, sid])
    Sigmaj[dif, dif] = SigRS + Sigmaj[dif, sid] %*% solve(Sigmaj[sid, sid]) %*% Sigmaj[sid, dif]
  }

  # Next, completion operation for sampled variance matrix
  H = c1
  for(i in 2:numberofcliques){
    dif = sort(setdiff(Cliques[[i]], H))
    sid = Separators[[i]]
    h = sort(setdiff(H, sid))
    Sigmaj[dif, h] = Sigmaj[dif, sid] %*% solve(Sigmaj[sid, sid]) %*% Sigmaj[sid, h]
    Sigmaj[h, dif] = t(Sigmaj[dif, h])
    H = sort(union(H, Cliques[[i]])) # probably no need to sort, just playing safe
  }
  Sigma = Sigmaj

  # Next, computing the corresponding sampled precision matrix
  Caux = Saux = array(0, c(p, p, numberofcliques))
  cid = Cliques[[1]]
  Caux[cid, cid, 1] = solve(Sigmaj[cid, cid])
  for(i in 2:numberofcliques){
    cid = Cliques[[i]]
    Caux[cid, cid, i] = solve(Sigmaj[cid, cid])
    sid = Separators[[i]]
    Saux[sid, sid, i] = solve(Sigmaj[sid, sid])
  }
  Omega = rowSums(Caux, dims = 2) - rowSums(Saux, dims = 2)

  return(list(Sigma = Sigma, Omega = Omega))
}


# Input:  an adjacency  matrix A of a decomposable graph G
# Output: cell array G containing the cliques and separators of G
#         nodeIDs and nodenames are optional inputs

makedecompgraph = function(Adj){
  # first check if "Adj" is a matrix object
  if(is.matrix(Adj)==FALSE) { stop("the input must be a matrix object!") }
  # check if "Adj" is a square matrix
  if(dim(Adj)[1]!=dim(Adj)[2]) { stop("the input must be a square matrix") }
  # check if "Adj" is symmetric
  if(isSymmetric.matrix(Adj)==FALSE) { stop("the input must be a symmetric matrix") }

  p = nrow(Adj)
  Adj[Adj!=0] = 1 # set all non-zero entries of Adj to 1
  Adj = Adj - diag(diag(Adj)) + diag(p) # set all diagonal entries of Adj to be 1
  Order = 1:p
  i = 1
  Adj0 = Adj
  while(i<p){
    nn = apply(Adj[1:i, (i+1):p, drop=F], 2, sum)
    b = which.max(nn)
    Order[c(i+1, b+i)] = Order[c(b+i, i+1)]
    i = i + 1
    Adj = Adj0
    Adj = Adj[Order, Order]
  }

  numberofcliques = 1
  Cliques = list(1)
  i = 2
  while(i<=p){
    if( sum(Adj[i,Cliques[[numberofcliques]]])==length(Cliques[[numberofcliques]]) ){
      Cliques[[numberofcliques]] = c(Cliques[[numberofcliques]], i)
    }
    else{
      numberofcliques = numberofcliques + 1
      Cliques[[numberofcliques]] = union(i, which(Adj[i, 1:i]==1))
    }
    i = i + 1
  }

  for(i in 1:numberofcliques){
    Cliques[[i]] = sort(Order[Cliques[[i]]])
  }

  UN = Cliques[[1]]
  Separators = list()
  if(numberofcliques==1){ return(list(C = Cliques, S = Separators)) }
  else{
    for(i in 2:numberofcliques){
      Separators[[i]] = sort(intersect(UN, Cliques[[i]]))
      UN = union(UN, Cliques[[i]])
    }
    return(list(C = Cliques, S = Separators))
  }
}





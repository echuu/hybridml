
source("examples/gwish/gwish_density.R")
library(BDgraph)
#### initialize graphs ---------------------------------------------------------



### create parameters for a single G_5 -----------------------------------------
p1 = 5
G_5 = matrix(c(1,1,0,1,1,
               1,1,1,0,0,
               0,1,1,1,1,
               1,0,1,1,1,
               1,0,1,1,1), p1, p1)

b = 100
n = 100
V_5 = n * diag(1, p1)

FREE_PARAMS_ALL = c(upper.tri(diag(1, p), diag = T) & G_5)
edgeInd = G_5[upper.tri(G_5, diag = TRUE)] %>% as.logical

## construct A matrix so that we can compute k_i
A = (upper.tri(diag(1, p1), diag = F) & G_5) + 0

k_i  = colSums(A) # see step 2, p. 329 of Atay
nu_i = rowSums(A) # see step 2, p. 329 of Atay
b_i = nu_i + k_i + 1
b_i

# S = matrix(0, p, p)
D = sum(edgeInd) # number of free parameters / dimension of parameter space

index_mat = matrix(0, p1, p1)
index_mat[upper.tri(index_mat, diag = T)][edgeInd] = 1:D
index_mat[upper.tri(index_mat, diag = T)]
t_ind = which(index_mat!=0,arr.ind = T)
t_ind

index_mat[lower.tri(index_mat)] = NA
vbar = which(index_mat==0,arr.ind = T) # non-free elements
n_nonfree = nrow(vbar)

params_G5 = list(G = G_5, P = P, p = p1, D = D, edgeInd = edgeInd,
                 b = b, nu_i = nu_i, b_i = b_i,
                 t_ind = t_ind, n_nonfree = n_nonfree, vbar = vbar)


################################################################################

### create parameters for stacked G_5's so we can sample, and compute psi ------

n_G5 = 3 # number of G_5 graphs we want to stack
G = diag(1, n_G5) %x% G_5
p = ncol(G)
V = n * diag(1, p)

## try computing the normalizing constant of G_9 first as sanity check
FREE_PARAMS_ALL = c(upper.tri(diag(1, p), diag = T) & G)

edgeInd = G[upper.tri(G, diag = TRUE)] %>% as.logical

## construct A matrix so that we can compute k_i
A = (upper.tri(diag(1, p), diag = F) & G) + 0

k_i  = colSums(A) # see step 2, p. 329 of Atay
nu_i = rowSums(A) # see step 2, p. 329 of Atay
b_i = nu_i + k_i + 1
b_i

set.seed(1)
Omega_G = rgwish(1, G, b, V) # generate the true precision matrix
P = chol(solve(V)) # upper cholesky factor; D^(-1) = TT'  in Atay paper

# params = list(G = G, P = P, p = p, edgeInd = edgeInd,
#               b = b, nu_i = nu_i, b_i = b_i)
N = 0
S = matrix(0, p, p)
D = sum(edgeInd) # number of free parameters / dimension of parameter space

index_mat = matrix(0, p, p)
index_mat[upper.tri(index_mat, diag = T)][edgeInd] = 1:D
index_mat[upper.tri(index_mat, diag = T)]
t_ind = which(index_mat!=0,arr.ind = T)
t_ind

index_mat[lower.tri(index_mat)] = NA
vbar = which(index_mat==0,arr.ind = T) # non-free elements
n_nonfree = nrow(vbar)

params = list(G = G, P = P, p = p, D = D, edgeInd = edgeInd,
              b = b, nu_i = nu_i, b_i = b_i,
              t_ind = t_ind, n_nonfree = n_nonfree, vbar = vbar)


samps = samplegw(J, G, b, N, V, S, solve(P), FREE_PARAMS_ALL)
u_samps = samps$Psi_free %>% data.frame
u_df = preprocess(u_samps, D, params)     # J x (D_u + 1)
u_df %>% head





source("C:/Users/ericc/Documents/hybridml/examples/gwish/gwish_density.R")
library(BDgraph)
library(dplyr)
library(hybridml)


### test code on diagonal scale matrix
#### initialize graphs ---------------------------------------------------------
### create parameters for a single G_5 -----------------------------------------
p = 5
G = matrix(c(1,1,0,1,1,
               1,1,1,0,0,
               0,1,1,1,1,
               1,0,1,1,1,
               1,0,1,1,1), p, p)
b = 5
n = 10
V = n * diag(1, p)
V = rgwish(1, G, b, diag(p))

P = chol(solve(V)) # upper cholesky factor; D^(-1) = TT'  in Atay paper

FREE_PARAMS_ALL = c(upper.tri(diag(1, p), diag = T) & G)
edgeInd = G[upper.tri(G, diag = TRUE)] %>% as.logical

## construct A matrix so that we can compute k_i
A = (upper.tri(diag(1, p), diag = F) & G) + 0

k_i  = colSums(A) # see step 2, p. 329 of Atay
nu_i = rowSums(A) # see step 2, p. 329 of Atay
b_i = nu_i + k_i + 1
b_i

# S = matrix(0, p, p)
D = sum(edgeInd) # number of free parameters / dimension of parameter space

index_mat = matrix(0, p, p)
index_mat[upper.tri(index_mat, diag = T)][edgeInd] = 1:D
index_mat[upper.tri(index_mat, diag = T)]
t_ind = which(index_mat!=0,arr.ind = T)
t_ind

index_mat[lower.tri(index_mat)] = NA
vbar = which(index_mat==0,arr.ind = T) # non-free elements
n_nonfree = nrow(vbar)
S = matrix(0, p, p)
params = list(G = G, P = P, p = p, D = D, edgeInd = edgeInd,
              b = b, nu_i = nu_i, b_i = b_i,
              t_ind = t_ind, n_nonfree = n_nonfree, vbar = vbar,
              FREE_PARAMS_ALL = FREE_PARAMS_ALL,
              n_graphs = 1)

## sample from posterior

set.seed(1)
J = 2000
samps = samplegw(J, G, b, 0, V, S, solve(P), FREE_PARAMS_ALL)
u_samps = samps$Psi_free %>% data.frame
# u_df = preprocess(u_samps, D, params)     # J x (D_u + 1)
# u_df = hybridml::gwish_preprocess(u_samps, D, params) # J x (D_u + 1)


Rcpp::sourceCpp("C:/Users/ericc/Documents/hybridml/examples/gwish/gwish.cpp")

# grad_gwish() <- dspi_ij <- dpsi         (updated version)
# grad_cpp() <- dpsi_cpp <- dpsi_rsij     (old version)

# hess_cpp() <- d2psi_ii, dpsi_rsij, d2 <- d2_rs         (old version)
# hess_gwish() <- d2psi_ijkl <- d2psi(r, s, i, j, k, l)  (updated version)

u = u_df[15,1:D] %>% unlist %>% unname
u_mat = create_psi_mat_cpp(u, params)


data.frame(grad_diag = grad_cpp(u, params),
           grad_numer = pracma::grad(f = psi, u, params = params),
           grad_general = grad_gwish(u, params),
           t_ind - 1)

library(microbenchmark)
microbenchmark(grad_numer = pracma::grad(f = psi, u, params = params),
               grad_general = grad_gwish(u, params),
               times = 10)

vbar - 1 # nonfree elements --> these are the ones that will have recursive calls
## test first row gradient calculation
rr = 0
ss = 2
ii = 0
jj = 0

dpsi(rr, ss, ii, jj, u_mat, params)
dpsi_ij(1, 1, u_mat, params)

## compare dpsi() output to dpsi_rsij() output
vbar - 1
rr = 1
ss = 4
ii = 1
jj = 1
dpsi(rr, ss, ii, jj, u_mat, params)
dpsi_rsij(rr, ss, ii, jj, u_mat, G)

data.frame(closed = grad_cpp(u, params),
           numerical = pracma::grad(f = psi, u, params = params))

### testing analytical vs. numerical hessian diagonal elements
data.frame(closed = diag(hess_cpp(u, params)),
           numerical = diag(pracma::hessian(psi, u, params = params)),
           closed_test = diag(hess_cpp_test(u, params)),
           # closed_general = diag(grad_gwish(u, params)), # this line crashes
           t_ind)

h_closed = hess_cpp(u, params)
h_numer = pracma::hessian(psi, u, params = params)
h_closed_update = hess_cpp_test(u, params)

all.equal(h_closed[upper.tri(h_closed)],
          h_numer[upper.tri(h_numer)])

all.equal(h_numer[upper.tri(h_closed)],
          h_closed_update[upper.tri(h_numer)])

all.equal(h_closed[upper.tri(h_closed)],
          h_closed_update[upper.tri(h_numer)])


### test code on non-diagonal scale matrix -------------------------------------

set.seed(2021)
p = 5; n = 300
G = matrix(0, p, p)
b = 3
G[1, p] = 1
for(i in 1:(p-1)){
  G[i, i+1] = 1
}

G = G + t(G); diag(G) = 1
Sigma = rgwish(1, G, b, diag(p))
Y = matrix(rnorm(n*p,), n, p)
Y = Y %*% t(chol(Sigma))
gnorm(G, n+b, diag(p)+t(Y)%*%Y, 1000)

















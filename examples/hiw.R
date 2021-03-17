

source("examples/hiw_helper.R")

testG  = matrix(c(1,1,0,0,0,
                  1,1,1,1,1,
                  0,1,1,1,1,
                  0,1,1,1,1,
                  0,1,1,1,1), 5, 5)

D = nrow(testG)
b = 3          # prior degrees of freedom
V = diag(1, D) # prior scale matrix

D_0 = 0.5 * D * (D + 1) # num entries on diagonal and upper diagonal
J = 1000

# logical vector determining existence of edges between vertices
edgeInd = testG[upper.tri(testG, diag = TRUE)] %>% as.logical
upperInd = testG[upper.tri(testG)] %>% as.logical
D_u = sum(edgeInd)

# Specify true value of Sigma
set.seed(1)
true_params = HIWsim(testG, b, V)
Sigma_G = true_params$Sigma
Omega_G = true_params$Omega # precision matrix -- this is the one we work with

# chol(Omega_G)

# Generate data Y based on Sigma_G
N = 100
Y = matrix(0, N, D)
for (i in 1:N) {
  Y[i, ] = t(t(chol(Sigma_G)) %*% rnorm(D, 0, 1)) # (500 x D)
}

S = t(Y) %*% Y

params = list(N = N, D = D, D_0 = D_0, testG = testG, edgeInd = edgeInd,
              upperInd = upperInd, S = S, V = V, b = b)

J = 5000
postIW = sampleHIW(J, D_u, D_0, testG, b, N, V, S, edgeInd)
post_samps = postIW$post_samps                 # (J x D_u)

u_df = hybridml::preprocess(post_samps, D_u, params)     # J x (D_u + 1)

(LIL = logmarginal(Y, testG, b, V, S))

lambda = function(u, params) { pracma::grad(psi, u, params = params) }
hess   = function(u, params) { pracma::hessian(psi, u, params = params) }

hybridml::hybml_const(u_df)$zhat
hybridml::hybml(u_df, params, grad = lambda, hess = hess)

(LIL = logmarginal(Y, testG, b, V, S))

- 0.5 * D * N * log(2 * pi) + BDgraph::gnorm(testG, b + N, V + S, iter = 1000) -
  BDgraph::gnorm(testG, b, V, iter = 1000)


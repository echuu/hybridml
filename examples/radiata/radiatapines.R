

source("examples/radiata/ModelEvidence.R")
source("examples/radiata/radiata_helper.R")
Rcpp::sourceCpp("examples/radiata/radiata.cpp")
RadiataPine = read.table("examples/radiata/RadiataPine.txt",
                         sep = " ", header = TRUE)

## data
y  = RadiataPine$y
n  = length(y)
X1 = cbind(rep(1,n),RadiataPine$x - mean(RadiataPine$x))
X2 = cbind(rep(1,n),RadiataPine$z - mean(RadiataPine$z))
d  = ncol(X1)

# normal prior on (alpha, beta)  ~ N (mu_0, (tau*Lamba0)^(-1))       -- model 1
# normal prior on (gamma, delta) ~ N (mu_0, (lambda * Lambda0)^(-1)) -- model 2
mu0 = c(3000,185)
Lambda0 = diag(1,2)
Lambda0[1,1] = 0.06
Lambda0[2,2] = 6.00

# gamma prior on
a_0 = 3
b_0 = 2*300^2

m1 = radiata_model(y, X1, mu0, Lambda0, a_0, b_0) # initialize model 1 object
m2 = radiata_model(y, X2, mu0, Lambda0, a_0, b_0) # initialize model 2 object

m1$logz # true log marginal likelihood under model 1
m2$logz # true log marginal likelihood under model 2

n.its   = 40000 # num iterations to run MCMC
burn.in = 10000 # num burn-in

m1_samps = m1$gibbs_radiata(n.its, burn.in, m1$fix)
m2_samps = m2$gibbs_radiata(n.its, burn.in, m2$fix)

u_df_m1 = hybridml::preprocess(m1_samps, d + 1, m1$params)
u_df_m2 = hybridml::preprocess(m2_samps, d + 1, m2$params)

## hybrid approximation (EP)
out1 = hybridml::hybml(u_df_m1, m1$params, grad = grad, hess = hess)
out2 = hybridml::hybml(u_df_m2, m2$params, grad = grad, hess = hess)

out1$logz
out2$logz

## bridge sampler estimate
log_density = function(u, data) {
  -psi(u, data)
}

u_samp = as.matrix(m1_samps)
colnames(u_samp) = names(u_df_m2)[1:(d+1)]
lb = c(rep(-Inf, d), 0)
ub = rep(Inf, d + 1)
names(lb) <- names(ub) <- colnames(u_samp)

bridge_result = bridgesampling::bridge_sampler(samples = u_samp,
                                               log_posterior = log_density,
                                               data = m1$params,
                                               lb = lb, ub = ub,
                                               silent = TRUE)
bridge_result$logml

m1$logz # true log marginal likelihood under model 1
m2$logz # true log marginal likelihood under model 2

exp(m2$logz - m1$logz)
exp(out2$logz - out1$logz)





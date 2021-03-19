
setwd("~/hybridml")

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
b_0 = 2 * 300^2

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
out1 = hybridml::hybml(u_df_m1, m1$params, grad = grad, hess = hess, u_0 = m1$theta_star())
out2 = hybridml::hybml(u_df_m2, m2$params, grad = grad, hess = hess, u_0 = m2$theta_star())

out1_map = hybridml::hybml(u_df_m1, m1$params, grad = grad, hess = hess)
out2_map = hybridml::hybml(u_df_m2, m2$params, grad = grad, hess = hess)

out1$logz
abs(m1$logz - out1$logz)

out2$logz
abs(m2$logz - out2$logz)

m1$logz # true log marginal likelihood under model 1
m2$logz # true log marginal likelihood under model 2

## true bayes factor
exp(m2$logz - m1$logz)

# hybrid bayes factor
exp(out2$logz - out1$logz)

# hybrid map bayes factor
exp(out2_map$logz - out1_map$logz)


n_reps = 20
lil_m1_map = numeric(n_reps)
lil_m1_opt = numeric(n_reps)
lil_m2_map = numeric(n_reps)
lil_m2_opt = numeric(n_reps)
bf_map     = numeric(n_reps)
bf_opt     = numeric(n_reps)


hybml(u_df_m1, m1$params, grad = grad, hess = hess, u_0 = m1$theta_star())

for (i in 1:n_reps) {

  m1_samps = m1$gibbs_radiata(n.its, burn.in, m1$fix)
  m2_samps = m2$gibbs_radiata(n.its, burn.in, m2$fix)

  u_df_m1 = hybridml::preprocess(m1_samps, d + 1, m1$params)
  u_df_m2 = hybridml::preprocess(m2_samps, d + 1, m2$params)

  ## hybrid approximation (EP)
  out1 = hybridml::hybml(u_df_m1, m1$params, grad = grad, hess = hess, u_0 = m1$theta_star())
  out2 = hybridml::hybml(u_df_m2, m2$params, grad = grad, hess = hess, u_0 = m2$theta_star())

  out1_map = hybridml::hybml(u_df_m1, m1$params, grad = grad, hess = hess)
  out2_map = hybridml::hybml(u_df_m2, m2$params, grad = grad, hess = hess)

  lil_m1_map[i] = out1_map$logz
  lil_m1_opt[i] = out1$logz
  lil_m2_map[i] = out2_map$logz
  lil_m2_opt[i] = out2$logz
  bf_map[i]     = exp(out2$logz - out1$logz)
  bf_opt[i]     = exp(out2_map$logz - out1_map$logz)

  print(paste("bf_map = ",   round(bf_map[i], 3),
              ", bf_opt = ", round(bf_opt[i], 3), sep = ''))

}

bf_21 = exp(m2$logz - m1$logz)
approx_df = data.frame(bf_map = bf_map, bf_opt = bf_opt, bf_21 = bf_21)
approx_df_long = reshape2::melt(approx_df, id.vars = 'bf_21')
x11()
ggplot(approx_df_long, aes(x = variable, y = value)) + geom_boxplot() +
  geom_hline(yintercept = bf_21, col = 'red', size = 1, linetype = 'dashed')

m1_df = data.frame(m1_map = lil_m1_map, m1_opt = lil_m1_opt, lil =  m1$logz)
m1_long = reshape2::melt(m1_df, id.vars = 'lil')
ggplot(m1_long, aes(x = variable, y = value)) + geom_boxplot() +
  geom_hline(yintercept = m1$logz, col = 'red', size = 1, linetype = 'dashed')

mean(bf_map)
mean(bf_opt)
bf_21




# ------------------------------------------------------------------------------

## bridge sampler estimate
log_density = function(u, data) {
  -psi(u, data)
}

u1_samp = as.matrix(m1_samps)
colnames(u1_samp) = names(u_df_m1)[1:(d+1)]
lb = c(rep(-Inf, d), 0)
ub = rep(Inf, d + 1)
names(lb) <- names(ub) <- colnames(u1_samp)

bridge_result1 = bridgesampling::bridge_sampler(samples = u1_samp,
                                                log_posterior = log_density,
                                                data = m1$params,
                                                lb = lb, ub = ub,
                                                silent = TRUE)

u2_samp = as.matrix(m2_samps)
colnames(u2_samp) = names(u_df_m2)[1:(d+1)]
lb = c(rep(-Inf, d), 0)
ub = rep(Inf, d + 1)
names(lb) <- names(ub) <- colnames(u2_samp)

bridge_result2 = bridgesampling::bridge_sampler(samples = u2_samp,
                                                log_posterior = log_density,
                                                data = m2$params,
                                                lb = lb, ub = ub,
                                                silent = TRUE)


b1_logz = bridge_result1$logml
b2_logz = bridge_result2$logml

# bridge bayes factor
exp(b2_logz - b1_logz)






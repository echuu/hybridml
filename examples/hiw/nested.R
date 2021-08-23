


N = 1000

logw = log(1.0 - exp(-1.0/N))

Theta = matrix(nrow = d,ncol = N)

log.likelihoods = numeric(N)

postIW = sampleHIW(N, D_u, D_0, G, b, N, V, S, edgeInd)
Theta = postIW$post_samps # (J x D_u)
u_df = hybridml::preprocess(post_samps, D_u, params)     # J x (D_u + 1)

log.likelihoods = apply(Theta, 2, cov_loglik, params = params)

## can compare above to the R version
logZ = -.Machine$double.xmax

log.contribution = 0

nest = 1

breakwhile = FALSE

#same termination criterion as Chopin and Robert (2010)

while((breakwhile==FALSE)){






} ## end nested sampling loop


logw = -nest/N - log(N)

for(i in 1:N){

  logZnew = evidence.obj$PLUS(logZ,logw+log.likelihoods[i])

  logZ = logZnew

}


logZ ## this is the estimate




radiata_model = function(y, X, mu0, Lambda0, a_0, b_0) {

  radiata_model_obj = list() # what

  tX = t(X)
  d = ncol(X)
  tau0 = Lambda0
  alpha = 2 * a_0
  delta = 2 * b_0

  # convenient constants
  logPi = log(pi)
  log2Pi = log(2*pi)
  XTX = tX %*% X
  XTy = tX %*% y
  M = XTX + tau0
  cholM = chol(M)
  log.detM = 2 * sum(log(diag(cholM)))
  invM = chol2inv(cholM)
  choltau0 = chol(tau0)
  invtau0 = chol2inv(choltau0)
  log.dettau0 = 2 * (sum(log(diag(choltau0))))
  P = diag(1,n) - X %*% invM %*% tX
  beta0 = invM %*% (tX %*% y + tau0 %*% mu0)
  yTy = t(y) %*% y
  c0 = yTy + t(mu0) %*% (tau0 %*% mu0) - t(beta0) %*% M %*% beta0
  c1 = t(mu0) %*% tau0 %*% mu0

  ## true log marginal likelihood
  LIL = -0.5*n*logPi + 0.5*log.dettau0 - 0.5*log.detM + 0.5*alpha*log(delta) +
    lgamma((n+alpha)/2) - lgamma(alpha/2) -0.5*(n+alpha)*log(c0+delta)


  params = list(Q_0 = Lambda0, mu0 = mu0, alpha = alpha, delta = delta,
                d = d, n = n, M = M,
                X = X, y = y, Xty = XTy,
                tau0 = Lambda0, beta0 = beta0,
                ldtau0 = log.dettau0)

  # list needed to pass into the gibbs sampler
  fix = list(); fix$vars = rep(FALSE, d + 1); fix$values = numeric(d + 1);


  radiata_model_obj$gibbs_radiata = function(Its, BurnIn, fix, initial = NULL,
                                             return.log.posterior = FALSE,
                                             return.log.likelihood = FALSE) {


    # print(head(X))

    # do site by site updates for fair comparison between methods
    T = matrix(nrow = Its - BurnIn, ncol = d + 1)

    # inialize from prior
    if (is.null(initial)) {
      tau = rgamma(1, shape = alpha / 2, rate = delta / 2)
      beta = rnorm(d, mean = mu0, sd = sqrt(1 / (tau * diag(tau0))))
    } else{
      beta = initial[1:d]
      tau = intial[d+1]
    }

    sh = 0.5 * (n + d + alpha)

    sample.vars = which(fix$vars[1:d] == FALSE)

    fix.vars = which(fix$vars[1:d] == TRUE)
    if(length(fix.vars) > 0) beta[fix.vars] = fix$values[fix.vars]

    if(fix$vars[d + 1] == TRUE) tau = fix$values[d + 1]

    sample.tau = !fix$vars[d+1]

    for(ItNum in 1:Its){

      #visit each parameter in turn
      for(j in sample.vars){
        w = M[j,]%*%(beta-beta0) - M[j,j]*(beta[j]-beta0[j])
        mu = beta0[j] - w/M[j,j]
        sig = sqrt(1/(tau*M[j,j]))
        beta[j] = rnorm(1,mean=mu,sd=sig)
      }

      rt = 0.5 * (t(beta-beta0) %*% M %*% (beta - beta0) + c0 + delta)
      if(sample.tau) tau = rgamma(1,shape = sh,rate = rt)

      if(ItNum > BurnIn){
        T[ItNum - BurnIn,] = c(beta, tau)
      }
    }

    return(as.data.frame(T))
  }

  radiata_model_obj$logz   = LIL
  radiata_model_obj$params = params
  radiata_model_obj$fix    = fix

  return(radiata_model_obj)
}

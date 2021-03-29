

epmgp_stable = function(m, K, b, lb, ub, EPS_CONVERGE = 1e-5) {

    D = length(lb)

    nu_site = numeric(D)
    tau_site = rep(0.5, D)

    ## initialize q(x)
    Sigma = K ## K is the covariance matrix that is passed in
    mu = (lb + ub) / 2

    Kinvm = solve(K) %*% m  ## don't need -- already computed this -> see below!
    Kinvm = b_k

    logZ = Inf             # store log normalizing constant
    mu_last = rep(-Inf, D) # store previous value of mu
    CONVERGED = FALSE      # monitor convergence of main algorithm loop
    k = 1

    ## initialize algorithm related quantities:
    tau_cavity = numeric(D)
    nu_cavity  = numeric(D)
    L          = matrix(0, D, D)
    logz_hat   = 0                # TODO: check dim
    sigma_hat  = numeric(D)       # TODO: check dim
    mu_hat     = numeric(D)       # TODO: check dim

    while (!CONVERGED) {

      # make cavity distribution
      tau_cavity = 1 / diag(Sigma) - tau_site
      nu_cavity  = mu / diag(Sigma) - nu_site

      if (any(tau_cavity < 0)) {
        print("sign error in tau_cavity -- flipping sign")
        tau_cavity = abs(tau_cavity)
        ## if this doesn't work, can look into updating tau_site only if
        ## it doesn't result in a negative estimate for tau_cavity
        ## although this looks like it results in really large values..
      }

      # moments = epmgp::trunc_norm_moments(lb, ub, nu_cavity / tau_cavity, 1 / tau_cavity)
      moments = tnorm_moments(lb, ub, nu_cavity / tau_cavity, 1 / tau_cavity)

      sigma_hat_out = abs(c(moments$sigma_hat))
      logz_hat_out  = c(moments$logz_hat)
      mu_hat_out    = c(moments$mu_hat)
      logz_hat      = c(logz_hat_out)
      sigma_hat     = c(sigma_hat_out)
      mu_hat        = c(mu_hat_out)

      # update site parameters
      delta_tau_site  = 1 / sigma_hat - tau_cavity - tau_site
      tau_site        = tau_site + delta_tau_site
      nu_site         = (mu_hat / sigma_hat) - nu_cavity

      # tau_site

      # enforce non-negativity of tau_site
      if (any(tau_site < 0)) {
        print("tau_site estimate < 0, setting to 0")
        tau_site[tau_site < 0] = 0
        # tau_site[tau_site > -1e-8] = 0
      }

      # update q(x) Sigma, mu
      S_site_half = diag(sqrt(tau_site))
      L = chol(diag(1, D) + S_site_half %*% K %*% S_site_half)
      V = solve(t(L), S_site_half %*% K)
      Sigma = K - t(V) %*% V
      mu = Sigma %*% (nu_site + Kinvm)

      if (pracma::Norm(mu_last - mu, 2) < EPS_CONVERGE) {
        # print(k)
        CONVERGED = TRUE
      } else {
        mu_last = mu
      }
      k = k + 1

    }

    ## compute logZ
    if (logZ != -Inf) {
      lZ1 = 0.5 * sum(log(1 + tau_site / tau_cavity)) - sum(log(diag(L)))
      lZ2 =  0.5 * t(nu_site - tau_site * m) %*%
        (Sigma - diag(1 / (tau_cavity + tau_site))) %*%
        (nu_site - tau_site * m)
      lZ3 = 0.5 * t(nu_cavity) %*% solve(diag(tau_site) + diag(tau_cavity),
                                         tau_site * nu_cavity / tau_cavity -
                                           2 * nu_site)
      lZ4 = -0.5 * t(tau_cavity * m) %*% solve(diag(tau_site) + diag(tau_cavity),
                                               tau_site * m - 2 * nu_site)
      logZ = lZ1 + lZ2 + lZ3 + lZ4 + sum(logz_hat)
    }

    result = list(logZ  = logZ,
                  mu    = mu,
                  Sigma = Sigma,
                  iter  = k)

    return(result)
}

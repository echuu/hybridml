
#include "gwish.h"
#include <cmath>
#define RCPP_ARMADILLO_RETURN_COLVEC_AS_VECTOR

typedef unsigned int u_int;

// [[Rcpp::depends(RcppArmadillo)]]
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS

/** utility functions **/
arma::vec chol2vec(arma::mat& M, u_int D);
arma::mat create_psi_mat_cpp(arma::vec u, Rcpp::List& params);


double psi_cpp_mat(arma::mat& psi_mat, Rcpp::List& params);
/** gradient related functions **/
arma::vec grad_cpp(arma::vec& u, Rcpp::List& params);
arma::vec grad_cpp_mat(arma::mat& psi_mat, Rcpp::List& params);
double dpsi_cpp(u_int i, u_int j, arma::mat& psi_mat, Rcpp::List& params);
double dpsi_rsij(u_int r, u_int s, u_int i, u_int j, 
	arma::mat& psi_mat, arma::mat& G);

/** hessian related functions **/
arma::mat hess_cpp(arma::vec& u, Rcpp::List& params);
arma::mat hess_cpp_mat(arma::mat& psi_mat, Rcpp::List& params);
double d2(u_int i, u_int j, u_int k, u_int l, 
	arma::mat& psi_mat, Rcpp::List& params);
double d2_rs(u_int r, u_int s, u_int i, u_int j, u_int k, u_int l,
	arma::mat& psi_mat, arma::mat& G);
double d2psi_ii(u_int r, u_int s, u_int i, arma::mat& psi_mat);

/** updated grad/hess functions for non-diagonal scale matrix **/
double xi(u_int i, u_int j, arma::mat& L);
arma::vec grad_gwish(arma::vec& u, Rcpp::List& params);
double dpsi_ij(u_int i, u_int j, arma::mat& psi_mat, Rcpp::List& params);
double dpsi(u_int r, u_int s, u_int i, u_int j, arma::mat& psi, Rcpp::List& params);


const double EPS_CONVERGE = 1e-8;


/** ------------------------------------------------------------------------ **/

// [[Rcpp::export]]
double erfcx (double x) {
  double a, d, e, m, p, q, r, s, t;
  
  a = fmax (x, 0.0 - x); // NaN preserving absolute value computation
  
  /* Compute q = (a-4)/(a+4) accurately. [0,INF) -> [-1,1] */
  m = a - 4.0;
  p = a + 4.0;
  r = 1.0 / p;
  q = m * r;
  t = fma (q + 1.0, -4.0, a); 
  e = fma (q, -a, t); 
  q = fma (r, e, q); 
  
  /* Approximate (1+2*a)*exp(a*a)*erfc(a) as p(q)+1 for q in [-1,1] */
  p =             0x1.edcad78fc8044p-31;  //  8.9820305531190140e-10
  p = fma (p, q,  0x1.b1548f14735d1p-30); //  1.5764464777959401e-09
  p = fma (p, q, -0x1.a1ad2e6c4a7a8p-27); // -1.2155985739342269e-08
  p = fma (p, q, -0x1.1985b48f08574p-26); // -1.6386753783877791e-08
  p = fma (p, q,  0x1.c6a8093ac4f83p-24); //  1.0585794011876720e-07
  p = fma (p, q,  0x1.31c2b2b44b731p-24); //  7.1190423171700940e-08
  p = fma (p, q, -0x1.b87373facb29fp-21); // -8.2040389712752056e-07
  p = fma (p, q,  0x1.3fef1358803b7p-22); //  2.9796165315625938e-07
  p = fma (p, q,  0x1.7eec072bb0be3p-18); //  5.7059822144459833e-06
  p = fma (p, q, -0x1.78a680a741c4ap-17); // -1.1225056665965572e-05
  p = fma (p, q, -0x1.9951f39295cf4p-16); // -2.4397380523258482e-05
  p = fma (p, q,  0x1.3be1255ce180bp-13); //  1.5062307184282616e-04
  p = fma (p, q, -0x1.a1df71176b791p-13); // -1.9925728768782324e-04
  p = fma (p, q, -0x1.8d4aaa0099bc8p-11); // -7.5777369791018515e-04
  p = fma (p, q,  0x1.49c673066c831p-8);  //  5.0319701025945277e-03
  p = fma (p, q, -0x1.0962386ea02b7p-6);  // -1.6197733983519948e-02
  p = fma (p, q,  0x1.3079edf465cc3p-5);  //  3.7167515521269866e-02
  p = fma (p, q, -0x1.0fb06dfedc4ccp-4);  // -6.6330365820039094e-02
  p = fma (p, q,  0x1.7fee004e266dfp-4);  //  9.3732834999538536e-02
  p = fma (p, q, -0x1.9ddb23c3e14d2p-4);  // -1.0103906603588378e-01
  p = fma (p, q,  0x1.16ecefcfa4865p-4);  //  6.8097054254651804e-02
  p = fma (p, q,  0x1.f7f5df66fc349p-7);  //  1.5379652102610957e-02
  p = fma (p, q, -0x1.1df1ad154a27fp-3);  // -1.3962111684056208e-01
  p = fma (p, q,  0x1.dd2c8b74febf6p-3);  //  2.3299511862555250e-01
  
  /* Divide (1+p) by (1+2*a) ==> exp(a*a)*erfc(a) */
  d = a + 0.5;
  r = 1.0 / d;
  r = r * 0.5;
  q = fma (p, r, r); // q = (p+1)/(1+2*a)
  t = q + q;
  e = (p - q) + fma (t, -a, 1.0); // residual: (p+1)-q*(1+2*a)
  r = fma (e, r, q);
  
  /* Handle argument of infinity */
  if (a > 0x1.fffffffffffffp1023) r = 0.0;
  
  /* Handle negative arguments: erfcx(x) = 2*exp(x*x) - erfcx(|x|) */
  if (x < 0.0) {
    s = x * x;
    d = fma (x, x, -s);
    e = exp (s);
    r = e - r;
    r = fma (e, d + d, r); 
    r = r + e;
    if (e > 0x1.fffffffffffffp1023) r = e; // avoid creating NaN
  }
  return r;
}

// [[Rcpp::export]]
Rcpp::List trunc_norm_moments(arma::vec lb_in, arma::vec ub_in, 
                              arma::vec mu_in, arma::vec sigma_in) {
  int d = lb_in.n_elem;
  arma::vec logz_hat_out(d);
  arma::vec z_hat_out(d);
  arma::vec mu_hat_out(d);
  arma::vec sigma_hat_out(d);
  
  for (int i = 0; i < d; i++) {
    
    double lb = lb_in(i);
    double ub = ub_in(i);
    double mu = mu_in(i);
    double sigma = sigma_in(i);
    
    double logz_hat_other_tail;
    double logz_hat;
    double mean_const;
    double var_const;
    
    // establish bounds
    double a = (lb - mu) / std::sqrt(2 * sigma);
    double b = (ub - mu) / std::sqrt(2 * sigma);

    /*
    Rcpp::Rcout << "a = " << a << std::endl;
    Rcpp::Rcout << "b = " << b << std::endl;
    
    if (abs(b) > abs(a)) {
      Rcpp::Rcout << "|b| > |a|" << std::endl;
    } else if (abs(a) == abs(b)) {
      Rcpp::Rcout << "|a| = |b|" << std::endl;
    }
    */

    // stable calculation
    
    // problem case 1
    if (std::isinf(a) && std::isinf(b)) {
      // check the sign
      if (copysign(1.0, a) == copysign(1.0, b)) {
        // integrating from inf to inf, should be 0
        logz_hat_out(i) = -arma::datum::inf;
        z_hat_out(i) = 0.0;
        mu_hat_out(i) = a;
        sigma_hat_out(i) = 0.0;
        continue;
      }
      else {
        logz_hat_out(i) = 0.0;
        z_hat_out(i) = 1.0;
        mu_hat_out(i) = mu;
        sigma_hat_out(i) = sigma;
        continue;
      }
    }
    
    // problem case 2
    else if (a > b) {
      // bounds are wrong, return 0 by convention
      logz_hat_out(i) = -arma::datum::inf;
      z_hat_out(i) = 0;
      mu_hat_out(i) = mu;
      sigma_hat_out(i) = 0;
      continue;
    }
    
    // typical case 1
    else if (a == -arma::datum::inf) {
      // integrating up to b
      if (b > 26.0) {
        // will be very close to 1
        logz_hat_other_tail = std::log(0.5) + std::log(erfcx(b)) - std::pow(b, 2);
        logz_hat = std::log1p(-std::exp(logz_hat_other_tail));
      }
      else {
        // b is small enough 
        logz_hat = std::log(0.5) + std::log(erfcx(-b)) - std::pow(b, 2);
      }
      
      mean_const = -2.0 / erfcx(-b);
      var_const = (-2.0 / erfcx(-b)) * (ub + mu);
    }
    
    // typical case 2
    else if (b == arma::datum::inf) {
      // Rcpp::Rcout << "handling unbounded upper constraint" << std::endl;
      // Rcpp::Rcout << "a: " << a << std::endl;
      // integrate from a to inf
      if (a < -26.0) {
        // will be very close to 1
        logz_hat_other_tail = std::log(0.5) + std::log(erfcx(-a)) - std::pow(a, 2);
        logz_hat = std::log1p(-std::exp(logz_hat_other_tail));
      }
      else {
        // should be stable
        logz_hat = std::log(0.5) + std::log(erfcx(a)) - std::pow(a, 2);
      }
      
      mean_const = 2.0 / erfcx(a);
      var_const = (2.0 / erfcx(a)) * (lb + mu);
    }
    
    // typical case 3
    else {
      // range from a to b, need stable exponent calculation
      double exp_a2b2 = std::exp(std::pow(a, 2) - std::pow(b, 2));
      double exp_b2a2 = std::exp(std::pow(b, 2) - std::pow(a, 2));
      
      // Rcpp::Rcout << "exp_a2b2: " << exp_a2b2 << std::endl;
      // Rcpp::Rcout << "exp_b2a2: " << exp_b2a2 << std::endl;


      if (copysign(1.0, a) == copysign(1.0, b)) {
        // exploit symmetry in problem to make calculations stable for erfcx
        double maxab = std::max(std::abs(a), std::abs(b));
        double minab = std::min(std::abs(a), std::abs(b));
        
        logz_hat = 
          std::log(0.5) - std::pow(minab, 2) + 
          std::log( std::abs( 
              std::exp( -(std::pow(maxab, 2) - std::pow(minab, 2))) * 
                erfcx(maxab) - 
                erfcx(minab) ) );
        
        double erfcx_a = erfcx(std::abs(a));
        double erfcx_b = erfcx(std::abs(b));

        // Rcpp::Rcout << "erfcx_a: " << erfcx_a << std::endl;
        // Rcpp::Rcout << "erfcx_b: " << erfcx_b << std::endl;

        mean_const = 2. * copysign(1.0, a) * (
          1 / (( erfcx_a - exp_a2b2 * erfcx_b )) -
            1 / (( exp_b2a2 * erfcx_a - erfcx_b ))
        );
        var_const = 2. * copysign(1.0, a) * (
          (lb + mu) / (erfcx_a - exp_a2b2 * erfcx_b) - 
            (ub + mu) / (exp_b2a2 * erfcx_a - erfcx_b)
        );
        // Rcpp::Rcout << "mean_const: " << mean_const << std::endl;
        // Rcpp::Rcout << "var_const: " << var_const << std::endl;
        // Rcpp::Rcout << "random: " << (exp_b2a2 * erfcx_a - erfcx_b) << std::endl;
      }
      
      else {
        // the signs are different, so b > a and b >= 0 and a <= 0
        if (std::abs(b) >= std::abs(a)) { 
          
          if (a >= -26.0) {
            // do things normally
            // Rcpp::Rcout << "first if" << std::endl;
            logz_hat = std::log(0.5) - std::pow(a, 2) + std::log( 
              erfcx(a) - std::exp(-(std::pow(b, 2) - std::pow(a, 2))) * erfcx(b)  
            );
            
            mean_const = 2 * (
              1 / (erfcx(a) - exp_a2b2 * erfcx(b)) -
                1 / (exp_b2a2 * erfcx(a) - erfcx(b))
            );
            var_const = 2 * (
              (lb + mu) / (erfcx(a) - exp_a2b2 * erfcx(b)) - 
                (ub + mu) / (exp_b2a2 * erfcx(a) - erfcx(b))
            );
          }
          
          else {
            // a is too small, so put in something close to 2 instead 

            logz_hat = std::log(0.5) - std::pow(b, 2) + std::log( 
              erfcx(-b) - std::exp(-(std::pow(a, 2) - std::pow(b, 2))) * erfcx(-a)
            );
            
            mean_const = 2 * (
              1 / (erfcx(a) - exp_a2b2 * erfcx(b)) -
                1 / (std::exp(std::pow(b, 2)) * 2 - erfcx(b))
            );
            var_const = 2 * (
              (lb + mu) / (erfcx(a) - exp_a2b2 * erfcx(b)) -
                (ub + mu) / (std::exp(std::pow(b, 2)) * 2 - erfcx(b))
            );
          }
        } // end first if()
        
        else {
          // abs(a) is bigger than abs(b), so reverse the calculation
          if (b <= 26.0) {
            // do things normally but mirrored across 0
            // Rcpp::Rcout << "here" << std::endl;
            logz_hat = std::log(0.5) - std::pow(b, 2) + std::log(
              erfcx(-b) - std::exp(-(std::pow(a, 2) - std::pow(b, 2))) * erfcx(-a)
            );
            
            mean_const = -2 * (
              1 / (erfcx(-a) - exp_a2b2 * erfcx(-b)) -
                1 / (exp_b2a2 * erfcx(-a) - erfcx(-b))
            );
            var_const = -2 * (
              (lb + mu) / (erfcx(-a) - exp_a2b2 * erfcx(-b)) -
                (ub + mu) / (exp_b2a2 * erfcx(-a) - erfcx(-b)) 
            );
          }
          else {
            // b is too big, put something close to 2 instead
            logz_hat = std::log(0.5) + std::log(
              2. - std::exp(-std::pow(a, 2)) * erfcx(-a) - std::exp(-std::pow(b, 2)) * erfcx(b)
            );
            
            mean_const = -2 * (
              1 / (erfcx(-a) - std::exp(std::pow(a, 2)) * 2) - 
                1 / (exp_b2a2 * erfcx(-a) - erfcx(-b))
            );
            var_const = -2 * (
              (lb + mu) / (erfcx(-a) - std::exp(std::pow(a, 2)) * 2) - 
                (ub + mu) / (exp_b2a2 * erfcx(-a) - erfcx(-b))
            );
          }
        }
      }
    }
    
    // Rcpp::Rcout << "Log z hat: " << logz_hat << std::endl;
    
    double z_hat = std::exp(logz_hat);
    double mu_hat = mu + mean_const * std::sqrt(sigma / (2 * arma::datum::pi));
    double sigma_hat = 
      sigma + var_const * sqrt(sigma / (2 * arma::datum::pi)) +
      std::pow(mu, 2) - std::pow(mu_hat, 2);      
    
    logz_hat_out(i) = logz_hat;
    z_hat_out(i) = z_hat;
    mu_hat_out(i) = mu_hat;
    sigma_hat_out(i) = sigma_hat;
  }
  
  Rcpp::List result = Rcpp::List::create(
    Rcpp::_["logz_hat"] = logz_hat_out,
    Rcpp::_["z_hat"] = z_hat_out,
    Rcpp::_["mu_hat"] = mu_hat_out,
    Rcpp::_["sigma_hat"] = sigma_hat_out
  );
  
  return result;
}


// [[Rcpp::export]]
double ep_logz(arma::vec m, arma::mat K, arma::vec lb, arma::vec ub) {
  
  arma::vec nu_site = arma::zeros(K.n_rows);
  arma::vec tau_site = arma::zeros(K.n_rows);
  
  // initialize q(x)
  arma::mat Sigma = K;
  arma::vec mu = (lb + ub) / 2;
  for (int i = 0; i < mu.n_elem; i++) {
    if (std::isinf(mu(i))) {
      mu(i) = copysign(1.0, mu(i)) * 100;
    }
  }
  
  arma::mat Kinvm = arma::solve(K, m);
  // Rcpp::Rcout << "Kinvm: " << Kinvm << std::endl;
  double logZ = arma::datum::inf;
  arma::vec mu_last = -arma::datum::inf * arma::ones(arma::size(mu));
  double converged = false;
  int k = 1;
  
  // algorithm loop
  arma::vec tau_cavity;
  arma::vec nu_cavity;
  arma::mat L;
  arma::vec logz_hat;
  arma::vec sigma_hat;
  arma::vec mu_hat;
  
  while (!converged) {
    
    // make cavity distribution
    tau_cavity = 1 / arma::diagvec(Sigma) - tau_site;
    nu_cavity = mu / arma::diagvec(Sigma) - nu_site;
    
    // compute moments using truncated normals
    Rcpp::List moments = trunc_norm_moments(lb, ub, nu_cavity / tau_cavity, 1 / tau_cavity);
    arma::vec sigma_hat_out = moments["sigma_hat"];
    arma::vec logz_hat_out = moments["logz_hat"];
    arma::vec mu_hat_out = moments["mu_hat"];
    logz_hat = logz_hat_out;
    sigma_hat = sigma_hat_out;
    mu_hat = mu_hat_out;

    // Rcpp::Rcout << "sigma_hat: " << sigma_hat_out << std::endl;
    
    // update the site parameters
    arma::vec delta_tau_site = (1 / sigma_hat) - tau_cavity - tau_site;
    tau_site += delta_tau_site;
    nu_site = (mu_hat / sigma_hat) - nu_cavity;
    
    // Rcpp::Rcout << "nu_site: " << nu_site.t() << std::endl;

    // enforce nonnegativity of tau_site
    if (arma::any(tau_site < 0)) {
      for (int i = 0; i < tau_site.n_elem; i++) {
        if (tau_site(i) > -1e-8) {
          tau_site(i) = 0.0;
        }
      }
    }
    
    // update q(x) Sigma and mu
    arma::mat S_site_half = arma::diagmat(arma::sqrt(tau_site));
    L = arma::chol(
      arma::eye(K.n_rows, K.n_cols) + S_site_half * K * S_site_half);
    arma::mat V = arma::solve(L.t(), S_site_half * K);
    Sigma = K - V.t() * V;
    mu = Sigma * (nu_site + Kinvm);
    
    // Rcpp::Rcout << "tau site: " << tau_site << std::endl;
    // Rcpp::Rcout << "tau cavity: " << tau_cavity << std::endl;
    // Rcpp::Rcout << "L: " << L << std::endl;
    
    // check convergence criteria
    if ((arma::norm(mu_last - mu)) < EPS_CONVERGE) {
      // Rcpp::Rcout << "converged: " << k << " iters" << std::endl;
      converged = true;
    } else {
      mu_last = mu;
    }
    k++;

    // if (k == 3) { break; }

  } // end while loop
  
  if (logZ != -arma::datum::inf) {
    double lZ1 = 0.5 * arma::sum(arma::log(1 + tau_site / tau_cavity)) - 
      arma::sum(arma::log(arma::diagvec(L)));

    double lZ2 = 0.5 * arma::as_scalar(
      (nu_site - tau_site % m).t() * 
      (Sigma - arma::diagmat(1 / (tau_cavity + tau_site))) * 
      (nu_site - tau_site % m)
    );

    double lZ3 = 0.5 * arma::as_scalar(
      nu_cavity.t() * 
        arma::solve(arma::diagmat(tau_site) + arma::diagmat(tau_cavity), 
                    tau_site % nu_cavity / tau_cavity - 2 * nu_site)
    );

    double lZ4 = -0.5 * arma::as_scalar(
      (tau_cavity % m).t() * 
        arma::solve(arma::diagmat(tau_site) + arma::diagmat(tau_cavity), 
                    tau_site % m - 2 * nu_site)
    );

    // Rcpp::Rcout << "tmp: " << (nu_site - tau_site % m).t() << std::endl;
    /*
    Rcpp::Rcout << "lz1: " << lZ1 << std::endl;
    Rcpp::Rcout << "lz2: " << lZ2 << std::endl;
    Rcpp::Rcout << "lz3: " << lZ3 << std::endl;
    Rcpp::Rcout << "lz4: " << lZ4 << std::endl;
    Rcpp::Rcout << "logzhat: " << logz_hat << std::endl;
    */
    logZ = lZ1 + lZ2 + lZ3 + lZ4 + arma::sum(logz_hat);



  }

  /*
  Rcpp::List result = Rcpp::List::create(
    Rcpp::_["logZ"] = logZ,
    Rcpp::_["mu"] = mu,
    Rcpp::_["Sigma"] = Sigma
  );
  */
  
  return logZ;
}

/** ------------------------------------------------------------------------ **/
// [[Rcpp::export]]
double lse(arma::vec arr, int count) 
{
   if(count > 0){
      double maxVal = arr(0);
      double sum = 0;

      for (int i = 1 ; i < count ; i++){
         if (arr(i) > maxVal){
            maxVal = arr(i);
         }
      }

      for (int i = 0; i < count ; i++){
         sum += exp(arr(i) - maxVal);
      }
      return log(sum) + maxVal;

   }
   else
   {
      return 0.0;
   }
}




// [[Rcpp::export]]
double approx_integral(u_int K, arma::mat& psi_df, arma::mat& bounds,
	Rcpp::List& params) {

	u_int D           = params["D"];    // dimension of parameter space

	arma::vec log_terms(K, arma::fill::zeros);
	arma::vec G_k(K, arma::fill::zeros);

	arma::mat H_k(D, D, arma::fill::zeros);
	arma::mat H_k_inv(D, D, arma::fill::zeros);
	arma::vec lambda_k(D, arma::fill::zeros);
	arma::vec b_k(D, arma::fill::zeros);
	arma::vec m_k(D, arma::fill::zeros);

	arma::vec lb(D, arma::fill::zeros);
	arma::vec ub(D, arma::fill::zeros);

	arma::vec u_k(D, arma::fill::zeros);

	for (u_int k = 0; k < K; k++) {

		// Rcpp::Rcout<< k << std::endl;

		u_k = arma::conv_to< arma::vec >::from(psi_df.submat(k, 0, k, D-1));

    arma::mat psi_mat = create_psi_mat_cpp(u_k, params);

		// Rcpp::Rcout<< u_k << std::endl;
		H_k = hess_cpp_mat(psi_mat, params);
		H_k_inv = inv(H_k);
		lambda_k = grad_cpp_mat(psi_mat, params);
		b_k = H_k * u_k - lambda_k;
		m_k = H_k_inv * b_k;

		/*
		if (k == (K-1)) {
			Rcpp::Rcout<< m_k << std::endl;
		}
		*/


		// TODO: extract the lower and upper bounds of the k-th partition
		for (u_int d = 0; d < D; d++) { 
			lb(d) = bounds.row(k)(2 * d);
			ub(d) = bounds.row(k)(2 * d + 1);
		}

		double val = 0;
		double sign;
		log_det(val, sign, H_k); 
		// Rcpp::Rcout << val << std::endl;

		// TODO: load the epmgp code into the same directory so that we can use
		// the EP code directly without having to go back into R env
		G_k(k) = ep_logz(m_k, H_k_inv, lb, ub);

		// Rcpp::Rcout<< psi_df(k, D) << std::endl;
		
		log_terms(k) = D / 2 * std::log(2 * M_PI) - 0.5 * val - psi_df(k, D) + 
			arma::dot(lambda_k, u_k) - 
			(0.5 * u_k.t() * H_k * u_k).eval()(0,0) + 
			(0.5 * m_k.t() * H_k * m_k).eval()(0,0) + G_k(k);

			// Rcpp::Rcout << (0.5 * u_k.t() * H_k * u_k).eval()(0,0) << std::endl;
			// Rcpp::Rcout << (0.5 + m_k.t() * H_k * m_k).eval()(0,0) << std::endl;

		// float x =  (m_k.t() * H_k * m_k).eval()(0,0);
		// Rcpp::Rcout << x << std::endl;

	} // end for() over k

	// TODO: find log-sum-exp function in arma
	return lse(log_terms, K);
} // end approx_integral() function



// [[Rcpp::export]]
arma::vec chol2vec(arma::mat& M, u_int D) {

	u_int k = 0; 
	u_int D_0 = D * (D + 1) / 2;
	arma::vec u(D_0);
	for (u_int c = 0; c < D; c++) {
		for (u_int r = 0; r <= c; r++) {
			u(k++) = M(r, c);	
		}
	}
	return u;
}


// [[Rcpp::export]]
arma::mat create_psi_mat_cpp(arma::vec u, Rcpp::List& params) {

	u_int p           = params["p"];    // dimension of the graph G
	u_int D           = params["D"];    // dimension of parameter space
	u_int b           = params["b"];    // degrees of freedom
	arma::vec nu_i    = params["nu_i"]; // see p. 329 of Atay (step 2)
	arma::vec b_i     = params["b_i"];  // see p. 329 of Atay (step 2)
	arma::mat P       = params["P"];    // upper cholesky factor of V_n
	arma::mat G       = params["G"];    // graph G 


	/* convert u into the matrix version with free-elements populated */
	// boolean vectorized version of G
	arma::vec G_bool  = params["FREE_PARAMS_ALL"];
	arma::uvec ids = find(G_bool); // extract indices of the free elements
	arma::vec u_prime(p * p, arma::fill::zeros);
	u_prime.elem(ids) = u;
	// Rcpp::Rcout << u_prime << std::endl;

	arma::mat u_mat(D, D, arma::fill::zeros);
	u_mat = reshape(u_prime, p, p);
	// Rcpp::Rcout << u_mat << std::endl;

	/* compute the non-free elmts in the matrix using the free elmts */

	// float test = sum(u_mat[i,i:(j-1)] * P[i:(j-1), j]);
	// arma::mat x1 = psi_mat.submat(0, 0, 0, 3);
	//Rcpp::Rcout << x1 << std::endl; 
	// Rcpp::Rcout << psi_mat.submat(0, 0, 3, 0) << std::endl;
	// Rcpp::Rcout << arma::dot(x1, psi_mat.submat(0, 0, 3, 0)) << std::endl;

	// u_mat should have all free elements in it
	double x0, tmp1, tmp2;
	for (u_int i = 0; i < p; i++) {
		for (u_int j = i; j < p; j++) {
			if (G(i,j) > 0) {
				continue; // free element, so these already have entries
			}
			if (i == 0) { // first row
				// TODO: double check this calculation
				u_mat(i,j) = -1/P(j,j) * arma::dot(u_mat.submat(i, i, i, j-1),
												   P.submat(i, j, j-1, j));
			} else {
				x0 = -1/P(j,j) * arma::dot(u_mat.submat(i, i, i, j-1),
										   P.submat(i, j, j-1, j));

				arma::vec tmp(i, arma::fill::zeros);
				for (u_int r = 0; r < i; r++) { 
					tmp1 = u_mat(r,i) + arma::dot(u_mat.submat(r, r, r, i-1),
											P.submat(r, i, i-1, i)) / P(i,i);
					tmp2 = u_mat(r,j) + arma::dot(u_mat.submat(r, r, r, j-1),
											P.submat(r, j, j-1, j)) / P(j,j);
					tmp(r) = tmp1 * tmp2;
				}

				u_mat(i,j) = x0 - 1 / u_mat(i,i) * arma::sum(tmp);
			}
		} // end inner for() over j
	} // end outer for() over i

	return u_mat;
} // end create_psi_mat_cpp() function



// [[Rcpp::export]]
double psi_cpp(arma::vec& u, Rcpp::List& params) {

	u_int p           = params["p"];    // dimension of the graph G
	u_int D           = params["D"];    // dimension of parameter space
	u_int b           = params["b"];    // degrees of freedom
	arma::vec nu_i    = params["nu_i"]; // see p. 329 of Atay (step 2)
	arma::vec b_i     = params["b_i"];  // see p. 329 of Atay (step 2)
	arma::mat P       = params["P"];    // upper cholesky factor of V_n
	arma::mat G       = params["G"];    // graph G 

	arma::mat psi_mat = create_psi_mat_cpp(u, params);

	double psi_u = p * std::log(2);
	for (u_int i = 0; i < p; i++) {
		psi_u += (b + b_i(i) - 1) * std::log(P(i, i)) + 
			(b + nu_i(i) - 1) * std::log(psi_mat(i, i)) -
			0.5 * std::pow(psi_mat(i,i), 2);
		for (u_int j = i + 1; j < p; j++) {
			psi_u += -0.5 * std::pow(psi_mat(i,j), 2);
		}
	}

	return -psi_u;
} // end of psi_cpp() function


// [[Rcpp::export]]
double psi_cpp_mat(arma::mat& psi_mat, Rcpp::List& params) {

	u_int p           = params["p"];    // dimension of the graph G
	u_int D           = params["D"];    // dimension of parameter space
	u_int b           = params["b"];    // degrees of freedom
	arma::vec nu_i    = params["nu_i"]; // see p. 329 of Atay (step 2)
	arma::vec b_i     = params["b_i"];  // see p. 329 of Atay (step 2)
	arma::mat P       = params["P"];    // upper cholesky factor of V_n
	arma::mat G       = params["G"];    // graph G 

	double psi_u = p * std::log(2);
	for (u_int i = 0; i < p; i++) {
		psi_u += (b + b_i(i) - 1) * std::log(P(i, i)) + 
			(b + nu_i(i) - 1) * std::log(psi_mat(i, i)) -
			0.5 * std::pow(psi_mat(i,i), 2);
		for (u_int j = i + 1; j < p; j++) {
			psi_u += -0.5 * std::pow(psi_mat(i,j), 2);
		}
	}
	return -psi_u;
}



// [[Rcpp::export]]
arma::vec grad_cpp(arma::vec& u, Rcpp::List& params) {


	arma::mat G       = params["G"]; // graph G represented as adjacency matrix
	u_int p           = params["p"]; // dimension of the graph G
	arma::vec edgeInd = params["edgeInd"];
	u_int D           = params["D"]; // dimension of parameter space
	// TODO: implement create_psi_mat() function later; for now, we pass it in
	// arma::mat psi_mat = vec2chol(u, p)

	arma::mat psi_mat = create_psi_mat_cpp(u, params);

	// initialize matrix that can store the gradient elements
	arma::mat gg(D, D, arma::fill::zeros);
	// populate the gradient matrix entry by entry
	for (u_int i = 0; i < p; i++) {
		for (u_int j = i; j < p; j++) {
			if (G(i,j) > 0) {
				gg(i,j) = dpsi_cpp(i, j, psi_mat, params);
			}
		}
	}
	// convert the matrix back into a vector and return only the entries
	// that have a corresponding edge in the graph 
	arma::vec grad_vec = chol2vec(gg, p);
	arma::uvec ids = find(edgeInd);

	return grad_vec.elem(ids);
}


// [[Rcpp::export]]
arma::vec grad_cpp_mat(arma::mat& psi_mat, Rcpp::List& params) {


  arma::mat G       = params["G"]; // graph G represented as adjacency matrix
  u_int p           = params["p"]; // dimension of the graph G
  arma::vec edgeInd = params["edgeInd"];
  u_int D           = params["D"]; // dimension of parameter space
  // TODO: implement create_psi_mat() function later; for now, we pass it in
  // arma::mat psi_mat = vec2chol(u, p)

  // arma::mat psi_mat = create_psi_mat_cpp(u, params);

  // initialize matrix that can store the gradient elements
  arma::mat gg(D, D, arma::fill::zeros);
  // populate the gradient matrix entry by entry
  for (u_int i = 0; i < p; i++) {
    for (u_int j = i; j < p; j++) {
      if (G(i,j) > 0) {
        gg(i,j) = dpsi_cpp(i, j, psi_mat, params);
      }
    }
  }
  // convert the matrix back into a vector and return only the entries
  // that have a corresponding edge in the graph 
  arma::vec grad_vec = chol2vec(gg, p);
  arma::uvec ids = find(edgeInd);

  return grad_vec.elem(ids);
}


// [[Rcpp::export]]
double dpsi_cpp(u_int i, u_int j, arma::mat& psi_mat, Rcpp::List& params) {

	arma::mat G     = params["G"];    // graph G represented as adjacency matrix
	u_int p         = params["p"];    // dimension of the graph G
	u_int b         = params["b"];    // degrees of freedom
	arma::vec nu_i  = params["nu_i"]; //

	double d_ij; // derivative of psi wrt psi_ij (summation over 
				// the derivatives wrt the free elements)

	if (G(i, j) == 0) {
		return 0;
	}
	if (i == j) { 

		d_ij = 0; 

		for (u_int r = 0; r < p; r++) {
			for (u_int s = r; s < p; s++) {
				if (G(r,s) == 0) {
					// if psi_rs == 0, no derivative calculation, skip
					if (psi_mat(r,s) == 0) {
						continue;
					}
					// otherwise: call the derivative function
					d_ij += psi_mat(r,s) * dpsi_rsij(r, s, i, j, psi_mat, G);
				} // end if for checking G[r,s]
			} // end loop over s
		} // end loop over r

		return d_ij - (b + nu_i(i) - 1) / psi_mat(i,i) + psi_mat(i,i);

	} else {

		d_ij = 0;

		for (u_int r = 0; r < p; r++) {
			for (u_int s = r; s < p; s++) {
				if (G(r,s) == 0) {
					if (psi_mat(r,s) == 0) {
						continue;
					}
					d_ij += psi_mat(r,s) * dpsi_rsij(r, s, i, j, psi_mat, G);
					// Rcpp::Rcout << "G[" << r+1 << ", " << s+1 << \
					// 		"] = " << G(r,s) << std::endl;
				}
			} // end loop over s
		} // end loop over r

	} // end if-else

    // only get to this return statement if we go through the else()
	return d_ij + psi_mat(i,j);

} // end dpsi() function



/* 
dpsi_rsij() : compute derivative d psi_rs / d psi_ij
*/
// [[Rcpp::export]]
double dpsi_rsij(u_int r, u_int s, u_int i, u_int j, 
	arma::mat& psi_mat, arma::mat& G) {

	// TODO: fill in the rest of the implementation
	if (G(r, s) > 0) { 
		if ((r == i) && (s == j)) { // d psi_{ij} / d psi_{ij} = 1
			return 1;
		} else { // d psi_{rs} / d psi_{ij} = 0, since psi_rs is free
			return 0;         
		}
	}

	if (i > r)                                { return 0; }
	if ((i == r) && (j > s))                  { return 0; }
	if ((i == r) && (j == s) && G(r, s) > 0)  { return 1; } // redundant check?
  	if ((i == r) && (j == s) && G(r, s) == 0) { return 0; } // d wrt to non-free

  	if ((r == 0) && (s > r)) { return 0; } // 1st row case -> simplified formula

  	if (r > 0) { 

  		// arma::vec tmp_sum(r - 1);
  		arma::vec tmp_sum(r);

  		for (u_int k = 0; k < r; k++) { 

  			if ((psi_mat(k,s) == 0) && (psi_mat(k,r) == 0)) {
  				tmp_sum(k) = 0;
  				continue;
  			} else {

  				if (psi_mat(k,s) == 0) {
  					tmp_sum(k) = -1/psi_mat(r,r) * 
  						dpsi_rsij(k, s, i, j, psi_mat, G) * psi_mat(k,r);
  				} else if (psi_mat(k,r) == 0) {
  					tmp_sum(k) = -1/psi_mat(r,r) * 
  						dpsi_rsij(k, r, i, j, psi_mat, G) * psi_mat(k,s);
  				} else {

  					if ((i == j) && (r == i) && (G(r,s) == 0)) {
  						tmp_sum(k) = 1/std::pow(psi_mat(r,r), 2) * 
  						psi_mat(k,r) * psi_mat(k,s) -
              				1/psi_mat(r,r) * 
              				(dpsi_rsij(k, r, i, j, psi_mat, G) * psi_mat(k, s) +
                  			 dpsi_rsij(k, s, i, j, psi_mat, G) * psi_mat(k, r));
			        } else {
			        	tmp_sum(k) = -1/psi_mat(r,r) * (
			              dpsi_rsij(k, r, i, j, psi_mat, G) * psi_mat(k, s) +
			              dpsi_rsij(k, s, i, j, psi_mat, G) * psi_mat(k,r));
			        }
  				}
  			} // end if-else
  		} // end for

  		// expression derived from Eq. (99)
  		return arma::sum(tmp_sum);

  	} else {
  		return -999;
  	}

} // end dpsi_rsij() function


// [[Rcpp::export]]
arma::mat hess_cpp_test(arma::vec& u, Rcpp::List& params) {

  u_int D           = params["D"];          // dimension of parameter space
  arma::mat G       = params["G"];          // graph G 
  u_int n_nonfree   = params["n_nonfree"];  // # nonfree elements
  arma::mat ind_mat = params["t_ind"];      // index of the free elements
  arma::mat vbar    = params["vbar"];       // index of nonfree elements
  u_int b           = params["b"];          // degrees of freedom
  arma::vec nu_i    = params["nu_i"];       //
  arma::mat psi_mat = create_psi_mat_cpp(u, params);

  arma::mat H(D, D, arma::fill::zeros);

  u_int d, i, j, a, r, c, rr, ss, k, l; // initialize various indices
  double tmp;

  // first populate the diagonal elements of the hessian matrix
  for (d = 0; d < D; d++) {
    // subtract one to account for 0-index in C++
    i = ind_mat(d, 0) - 1; // row loc of d-th free element
    j = ind_mat(d, 1) - 1; // col loc of d-th free element

    if (i == j) { // diagonal free elements 
      tmp = 0;
      for (a = 0; a < n_nonfree; a++) {

        rr = vbar(a, 0) - 1; // row index of non-free element
        ss = vbar(a, 1) - 1; // col index of non-free element

        tmp += std::pow(dpsi_rsij(rr, ss, i, j, psi_mat, G), 2) +
          psi_mat(rr,ss) * d2psi_ii(rr, ss, i, psi_mat);
      } // end of iteration over non-free elements

      H(d,d) = (b + nu_i(i) - 1) / std::pow(psi_mat(i,i), 2) + 1 + tmp;
    } else { 

      tmp = 0;
      for (a = 0; a < n_nonfree; a++) {
        rr = vbar(a, 0) - 1; // row index of non-free element
        ss = vbar(a, 1) - 1; // col index of non-free element

        // 11/5/21: previous implementation (line commented out below) 
        // did not account for the 2nd order derivative term
        // tmp += std::pow(dpsi_rsij(rr, ss, i, j, psi_mat, G), 2)
        tmp += std::pow(dpsi_rsij(rr, ss, i, j, psi_mat, G), 2) + 
          psi_mat(rr, ss) * d2_rs(rr, ss, i, j, i, j, psi_mat, G);
      }

      H(d,d) = 1 + tmp;
    } // end if-else
  } // end for() over d


  // populate the off-diagonal elements of H
  for (r = 0; r < (D-1); r++) { // index should be correct now

    i = ind_mat(r,0) - 1;
    j = ind_mat(r,1) - 1;

    for (c = r + 1; c < D; c++) { // index should be correct

      k = ind_mat(c, 0) - 1;
      l = ind_mat(c, 1) - 1;

      // H(r,c) = d2(i, j, k, l, psi_mat, params);
      // 11/5/21: testing to account for missing first order partial term
      H(r,c) = dpsi_rsij(i, j, k, l, psi_mat, G) + 
                d2(i, j, k, l, psi_mat, params);
      H(c,r) = H(r,c); // reflect calculation across diagonal
    } // end inner for() over cols
  } // end outer for() over rows
  return H;
}



// [[Rcpp::export]]
double d2psi_ii(u_int r, u_int s, u_int i, arma::mat& psi_mat) {

     double out = 0;
    // if (i == r)
    if (i == r) { 
      for (u_int m = 0; m < r; m++) { 
        out += psi_mat(m,r) * psi_mat(m,s);
      }


    } else {
      return out;
    }


    return(-2 / std::pow(psi_mat(r,r), 3) * out);
}


// [[Rcpp::export]]
arma::mat hess_cpp(arma::vec& u, Rcpp::List& params) {

	u_int D           = params["D"];          // dimension of parameter space
	arma::mat G       = params["G"];          // graph G 
	u_int n_nonfree   = params["n_nonfree"];  // # nonfree elements
	arma::mat ind_mat = params["t_ind"];      // index of the free elements
	arma::mat vbar    = params["vbar"];       // index of nonfree elements
	u_int b           = params["b"];          // degrees of freedom
	arma::vec nu_i    = params["nu_i"];       //
	arma::mat psi_mat = create_psi_mat_cpp(u, params);

	arma::mat H(D, D, arma::fill::zeros);

	u_int d, i, j, a, r, c, rr, ss, k, l; // initialize various indices
	double tmp;

	// first populate the diagonal elements of the hessian matrix
	for (d = 0; d < D; d++) {
		// subtract one to account for 0-index in C++
		i = ind_mat(d, 0) - 1; // row loc of d-th free element
		j = ind_mat(d, 1) - 1; // col loc of d-th free element

		if (i == j) {
			tmp = 0;
			for (a = 0; a < n_nonfree; a++) {

				rr = vbar(a, 0) - 1; // row index of non-free element
				ss = vbar(a, 1) - 1; // col index of non-free element

				tmp += std::pow(dpsi_rsij(rr, ss, i, j, psi_mat, G), 2);
			} // end of iteration over non-free elements

			H(d,d) = (b + nu_i(i) - 1) / std::pow(psi_mat(i,i), 2) + 1 + tmp;
		} else { 

			tmp = 0;
			for (a = 0; a < n_nonfree; a++) {
				rr = vbar(a, 0) - 1; // row index of non-free element
				ss = vbar(a, 1) - 1; // col index of non-free element

				tmp += std::pow(dpsi_rsij(rr, ss, i, j, psi_mat, G), 2);
			}

			H(d,d) = 1 + tmp;
		} // end if-else
	} // end for() over d


	// populate the off-diagonal elements of H
	for (r = 0; r < (D-1); r++) { // index should be correct now

		i = ind_mat(r,0) - 1;
		j = ind_mat(r,1) - 1;

		for (c = r + 1; c < D; c++) { // index should be correct

			k = ind_mat(c, 0) - 1;
			l = ind_mat(c, 1) - 1;

			H(r,c) = d2(i, j, k, l, psi_mat, params);
			H(c,r) = H(r,c); // reflect calculation across diagonal
		} // end inner for() over cols
	} // end outer for() over rows
	return H;
}

// [[Rcpp::export]]
arma::mat hess_cpp_mat(arma::mat& psi_mat, Rcpp::List& params) {

  u_int D           = params["D"];          // dimension of parameter space
  arma::mat G       = params["G"];          // graph G 
  u_int n_nonfree   = params["n_nonfree"];  // # nonfree elements
  arma::mat ind_mat = params["t_ind"];      // index of the free elements
  arma::mat vbar    = params["vbar"];       // index of nonfree elements
  u_int b           = params["b"];          // degrees of freedom
  arma::vec nu_i    = params["nu_i"];       //
  // arma::mat psi_mat = create_psi_mat_cpp(u, params);

  arma::mat H(D, D, arma::fill::zeros);

  u_int d, i, j, a, r, c, rr, ss, k, l; // initialize various indices
  double tmp;

  // first populate the diagonal elements of the hessian matrix
  for (d = 0; d < D; d++) {
    // subtract one to account for 0-index in C++
    i = ind_mat(d, 0) - 1; // row loc of d-th free element
    j = ind_mat(d, 1) - 1; // col loc of d-th free element

    if (i == j) {
      tmp = 0;
      for (a = 0; a < n_nonfree; a++) {

        rr = vbar(a, 0) - 1; // row index of non-free element
        ss = vbar(a, 1) - 1; // col index of non-free element

        tmp += std::pow(dpsi_rsij(rr, ss, i, j, psi_mat, G), 2);
      } // end of iteration over non-free elements

      H(d,d) = (b + nu_i(i) - 1) / std::pow(psi_mat(i,i), 2) + 1 + tmp;
    } else { 

      tmp = 0;
      for (a = 0; a < n_nonfree; a++) {
        rr = vbar(a, 0) - 1; // row index of non-free element
        ss = vbar(a, 1) - 1; // col index of non-free element

        tmp += std::pow(dpsi_rsij(rr, ss, i, j, psi_mat, G), 2);
      }

      H(d,d) = 1 + tmp;
    } // end if-else
  } // end for() over d


  // populate the off-diagonal elements of H
  for (r = 0; r < (D-1); r++) { // index should be correct now

    i = ind_mat(r,0) - 1;
    j = ind_mat(r,1) - 1;

    for (c = r + 1; c < D; c++) { // index should be correct

      k = ind_mat(c, 0) - 1;
      l = ind_mat(c, 1) - 1;

      H(r,c) = d2(i, j, k, l, psi_mat, params);
      H(c,r) = H(r,c); // reflect calculation across diagonal
    } // end inner for() over cols
  } // end outer for() over rows
  return H;
}


// [[Rcpp::export]]
double d2(u_int i, u_int j, u_int k, u_int l, 
	arma::mat& psi_mat, Rcpp::List& params) {

	u_int n_nonfree   = params["n_nonfree"];  // # nonfree elements
	arma::mat G       = params["G"];          // graph G 
	arma::mat vbar    = params["vbar"];       // index of nonfree elements

	arma::vec tmp(n_nonfree, arma::fill::zeros);

	u_int n, r, s;
	for (n = 0; n < n_nonfree; n++) {
		r = vbar(n, 0) - 1; // row index of nonfree element
		s = vbar(n, 1) - 1; // col index of nonfree element
		if (psi_mat(r,s) == 0) { // can be combined partially with step below
			tmp(n) = dpsi_rsij(r, s, k, l, psi_mat, G) * 
				dpsi_rsij(r, s, i, j, psi_mat, G);
		} else {
			// product rule: (d psi_rs / d psi_ij d psi_kl)^2 + 
      //               psi_rs * (d^2 psi_rs / d psi_ij ^ 2)
			tmp(n) = dpsi_rsij(r, s, k, l, psi_mat, G) * 
				dpsi_rsij(r, s, i, j, psi_mat, G) + 
				psi_mat(r,s) * d2_rs(r, s, i, j, k, l, psi_mat, G);
			// tmp(n) = d2_rs(r, s, i, j, k, l, psi_mat, G);
		} // end of if-else

	} // end for() over nonfree elements
	return arma::sum(tmp);
} // end d2() function


// [[Rcpp::export]]
double d2_rs(u_int r, u_int s, u_int i, u_int j, u_int k, u_int l,
	arma::mat& psi_mat, arma::mat& G) {

	/** assumption: we don't have i == k AND j == l, i.e., computing 2nd order
	    derivative for diagonal term in the hessian matrix, since we've already
	    done this previously in a loop dedicated to computing diagonal terms 
	    in the hessian matrix
	**/

	/*********** check for early break condtions to save computing ************/
	if (G(r, s) > 0) { return 0; } 			// free element
	if ((r < i) || (r < k)) { return 0; }   // row below
	// same row, col after
  // 11/4/21: condition below is too stringent, r doesn't need to equal BOTH i
  // AND k. Instead, we should check if either is equal and then provide
  // a check for their corresponding columns to see i
	// if ((r == i) && (r == k) && ((s < j) || (s < l))) { return 0; } 

  // 11/4/21: revised condition see notes to see explanation:
  // TODO: add in the explanation for this one later
  if (((r == i) && (s < j)) || ((r == k) && (s < l))) { return 0;}


	/********* general case: recursive call for 2nd order derivative **********/
	// see note below in loop for why vector is size r instead of size (r-1)
	arma::vec tmp(r, arma::fill::zeros); 
	u_int m;
	for (m = 0; m < r; m++) { 
		/* Note on the 0-indexing and upper limit of the loop:
		   Because r represents the *row number*, in the R code, r = 3 for row 3
		   so the loop will run 1:(3-1) -> 2 times, but in the C++ code, we do
		   not need to subtract 1 from the upper limit because row 3 corresponds
		   to r = 2, so the loop will run 0:1 -> 2 times, matching the number
		   of iterations run in the R code
		*/

		if ((psi_mat(m, s) == 0) && (psi_mat(m, r) == 0)) { 
			tmp(m) = - 1 / psi_mat(r, r) * (
				dpsi_rsij(m, r, i, j, psi_mat, G) * 
				dpsi_rsij(m, s, k, l, psi_mat, G) +
           		dpsi_rsij(m, r, k, l, psi_mat, G) * 
           		dpsi_rsij(m, s, i, j, psi_mat, G)
			);
		} else {

			if ((r == i) && (i == j)) { // case: d^2(psi_rs) / (dpsi_rr dpsi_kl)

				if (psi_mat(m, s) == 0) {
					tmp(m) = 1 / std::pow(psi_mat(r, r), 2) * 
						(psi_mat(m, r) * dpsi_rsij(m, s, k, l, psi_mat, G)) -
						1 / psi_mat(r, r) * (
						    d2_rs(m, s, i, j, k, l, psi_mat, G) * psi_mat(m,r) +
			               	dpsi_rsij(m, r, i, j, psi_mat, G) * 
			               	dpsi_rsij(m, s, k, l, psi_mat, G) +
			               	dpsi_rsij(m, r, k, l, psi_mat, G) * 
			               	dpsi_rsij(m, s, i, j, psi_mat, G)
               			);

				} else if (psi_mat(m, r) == 0) {
					tmp(m) = 1 / std::pow(psi_mat(r, r), 2) * 
						(dpsi_rsij(m, r, k, l, psi_mat, G) * psi_mat(m, s)) - 
						1 / psi_mat(r, r) * (
							d2_rs(m, r, i, j, k, l, psi_mat, G) * psi_mat(m,s) +
               				dpsi_rsij(m, r, i, j, psi_mat, G) * 
               				dpsi_rsij(m, s, k, l, psi_mat, G) +
               				dpsi_rsij(m, r, k, l, psi_mat, G) * 
               				dpsi_rsij(m, s, i, j, psi_mat, G)
						);

				} else {
					tmp(m) = 1 / std::pow(psi_mat(r, r), 2) * 
						(dpsi_rsij(m, r, k, l, psi_mat, G) * psi_mat(m, s) + 
						 psi_mat(m, r) * dpsi_rsij(m, s, k, l, psi_mat, G)) -
			            1 / psi_mat(r, r) * (
			            	d2_rs(m, r, i, j, k, l, psi_mat, G) * psi_mat(m,s) +
			               	d2_rs(m, s, i, j, k, l, psi_mat, G) * psi_mat(m,r) +
			               	dpsi_rsij(m, r, i, j, psi_mat, G) * 
			               	dpsi_rsij(m, s, k, l, psi_mat, G) +
			               	dpsi_rsij(m, r, k, l, psi_mat, G) * 
			               	dpsi_rsij(m, s, i, j, psi_mat, G)
			        	);
				} // 

			} else if ((r == k) && (k == l)) { 
				tmp(m) = d2_rs(r, s, k, l, i, j, psi_mat, G);
			} else { // case when r != i
				if (psi_mat(m, s) == 0) { 
					tmp(m) = -1 / psi_mat(r, r) * (
						dpsi_rsij(m,r,i,j,psi_mat,G) * 
						dpsi_rsij(m,s,k,l,psi_mat,G) +
               			dpsi_rsij(m,r,k,l,psi_mat,G) * 
               			dpsi_rsij(m,s,i,j,psi_mat,G) +
              			psi_mat(m,r) * d2_rs(m, s, i, j, k, l, psi_mat, G)
					);

				} else if (psi_mat(m, r) == 0) {
					tmp(m) = - 1 / psi_mat(r, r) * (
						dpsi_rsij(m, r, i, j, psi_mat, G) * 
						dpsi_rsij(m, s, k, l, psi_mat, G) +
               			dpsi_rsij(m, r, k, l, psi_mat, G) * 
               			dpsi_rsij(m, s, i, j, psi_mat, G) +
               			psi_mat(m, s) * d2_rs(m, r, i, j, k, l, psi_mat, G)
					);
				} else { 
					tmp(m) = - 1 / psi_mat(r, r) * (
						dpsi_rsij(m, r, i, j, psi_mat, G) * 
						dpsi_rsij(m, s, k, l, psi_mat, G) +
		               	dpsi_rsij(m, r, k, l, psi_mat, G) * 
		               	dpsi_rsij(m, s, i, j, psi_mat, G) +
		               	psi_mat(m, s) * d2_rs(m, r, i, j, k, l, psi_mat, G) +
		               	psi_mat(m, r) * d2_rs(m, s, i, j, k, l, psi_mat, G)
					);
				}

			}
		} // end of main if-else
	} // end for() over m

	return arma::sum(tmp);
} // end d2_rs() function


/**** update gradient, hessian to accommodate non-diagonal scale matrices *****/

// [[Rcpp::export]]
double xi(u_int i, u_int j, arma::mat& L) { 
    // L is UPPER TRIANGULAR cholesky factor of the INVERSE scale matrix, i.e,
    // D^(-1) = L'L
    return L(i, j) / L(j, j);
}


// [[Rcpp::export]]
arma::vec grad_gwish(arma::vec& u, Rcpp::List& params) {

    arma::mat G       = params["G"]; // graph G represented as adjacency matrix
    u_int p           = params["p"]; // dimension of the graph G
    arma::vec edgeInd = params["edgeInd"];
    u_int D           = params["D"]; // dimension of parameter space
    // TODO: implement create_psi_mat() function later; for now, we pass it in
    // arma::mat psi_mat = vec2chol(u, p)
    arma::mat psi_mat = create_psi_mat_cpp(u, params);

    // initialize matrix that can store the gradient elements
    arma::mat gg(D, D, arma::fill::zeros);
    // populate the gradient matrix entry by entry
    for (u_int i = 0; i < p; i++) {
        for (u_int j = i; j < p; j++) {
            if (G(i,j) > 0) {
                gg(i,j) = dpsi_ij(i, j, psi_mat, params);
            }
        }
    }
    // convert the matrix back into a vector and return only the entries
    // that have a corresponding edge in the graph 
    arma::vec grad_vec = chol2vec(gg, p);
    arma::uvec ids = find(edgeInd);

    return grad_vec.elem(ids);
} // end grad_gwish() function



// [[Rcpp::export]]
double dpsi_ij(u_int i, u_int j, arma::mat& psi_mat, Rcpp::List& params) {

    arma::mat G     = params["G"];    // graph G represented as adjacency matrix
    u_int p         = params["p"];    // dimension of the graph G
    u_int b         = params["b"];    // degrees of freedom
    arma::vec nu_i  = params["nu_i"]; //
    // arma::mat L     = params["P"];

    double d_ij; // derivative of psi wrt psi_ij (summation over 
                 // the derivatives wrt the free elements)
    if (G(i, j) == 0) {
        return 0;
    }
    if (i == j) { 
        d_ij = 0; 
        for (u_int r = 0; r < p; r++) {
            for (u_int s = r; s < p; s++) {
                if (G(r,s) == 0) {
                    // if psi_rs == 0, no derivative calculation, skip
                    if (psi_mat(r,s) == 0) {
                        continue;
                    }
                    // otherwise: call the derivative function
                    d_ij += psi_mat(r,s) * dpsi(r, s, i, j, psi_mat, params);
                } // end if for checking G[r,s]
            } // end loop over s
        } // end loop over r
        return d_ij - (b + nu_i(i) - 1) / psi_mat(i,i) + psi_mat(i,i);
    } else {
        d_ij = 0;
        for (u_int r = 0; r < p; r++) {
            for (u_int s = r; s < p; s++) {
                if (G(r,s) == 0) {
                    if (psi_mat(r,s) == 0) {
                        continue;
                    }
                    d_ij += psi_mat(r,s) * dpsi(r, s, i, j, psi_mat, params);
                    // Rcpp::Rcout << "G[" << r+1 << ", " << s+1 << \
                    //      "] = " << G(r,s) << std::endl;
                }
            } // end loop over s
        } // end loop over r
    } // end if-else
    // only get to this return statement if we go through the else()
    return d_ij + psi_mat(i,j);
} // end dpsi() function

// [[Rcpp::export]]
int test() { 

    double s0, s10, s11 = 0; 
    Rcpp::Rcout << "s0 = "  << s0 << std::endl;
    Rcpp::Rcout << "s10 = " << s10 << std::endl;
    Rcpp::Rcout << "s11 = " << s11 << std::endl;

    return 0;

}

// [[Rcpp::export]]
double dpsi(u_int r, u_int s, u_int i, u_int j, arma::mat& psi, Rcpp::List& params) {

    arma::mat G = params["G"];
    arma::mat L = params["P"];

    // Rcpp::Rcout << "beginning"  << std::endl;
    // Rcpp::Rcout << "G(" << r << ", " << s << ") = " << G(r,s) << std::endl;
    // Rcpp::Rcout << "r =  " << r << ", s =  " << s << ", i =  " << i << 
    //    ", j =  " << j << std::endl;

    if (G(r, s) > 0) { 
        if ((r == i) && (s == j)) { // d psi_{ij} / d psi_{ij} = 1
            return 1;
        } else { // d psi_{rs} / d psi_{ij} = 0, since psi_rs is free
            return 0;         
        }
    }
    if (i > r)                                { return 0; }
    if ((i == r) && (j > s))                  { return 0; }
    /* 11/5: fairly certain we don't need the following 2 checks */ 
    if ((i == r) && (j == s) && G(r, s) > 0)  { return 1; } // redundant check?
    if ((i == r) && (j == s) && G(r, s) == 0) { return 0; } // d wrt to non-free

    u_int k, l;
    double x;
    // 1st row case: simplified formulation of free elements
    if (r == 0) { 
        // don't need to check (s > r) because that case is already flagged by
        // the initial check in this function: if (G(r,s) > 0)
        x = 0;
        for (k = 0; k < s; k++) { 
            x += dpsi(r, k, i, j, psi, params) * xi(k, s, L); // correct G
        }
        return -x;
    } // end row 1 case

    bool DWRT_SAME_ROW_DIAG = ((i == j) && (r == i) && (G(r, s) == 0));
    // bool DWRT_SAME_ROW_DIAG = ((i == j) && (r == i));
    double s0 = 0, s10 = 0, s11 = 0; 
    double s12 = 0, s13 = 0, s110 = 0, s111 = 0, s120 = 0, s121 = 0; // i != j
    double out; // store the result of the calculation
    arma::vec s1(r); // store each term in the summation

    // TODO: THIS CHUNK IN PROGRESS
    if (DWRT_SAME_ROW_DIAG) { // dpsi_rs / dpsi_rr
        // Rcpp::Rcout << "HERE"  << std::endl;
        for (k = r; k < s; k++) { 
            s0 += xi(k, s, L) * dpsi(r, k, i, j, psi, params);
        }
        for (k = 0; k < r; k++) {
            for (l = k; l < s; l++) {
                s10 += psi(k,l) * xi(l, s, L);
            }
            for (l = k; l < r; l++) {
                s11 += psi(k,l) * xi(l, r, L);
            }
            s1(k) = psi(k,r) * psi(k,s) + 
                psi(k,r) * s10 + psi(k,s) * s11 + s10 * s11;
        } // end inner for()

        out = -s0 + 1 / std::pow(psi(r,r), 2) * arma::sum(s1);

    } else { // dpsi_rs / dpsi_ij
        // general case when derivative is wrt general (i, j), i > 0
        for (k = r; k < s; k++) {
            s0 += xi(k, s, L) * dpsi(r, k, i, j, psi, params);
        }
        /* calculation of s1, inner summation from 1:(r-1)
           s1 = s10 + s11 + s12 + s13, where each component is one of the four
           terms in the summation. s13 is consisted of 4 separate summations,
           for which we use 2 for loops to compute.
        */
        for (k = 0; k < r; k++) {

            // compute the intermediate summations: 
            for (l = k; l < s; l++) {
                s110 += psi(k, l) * xi(l, s, L);
                // note: s111 is just s110 w/ derivative wrt (ij) applied to 
                // the first term
                s111 += dpsi(k, l, i, j, psi, params) * xi(l, s, L); 
            }
            for (l = k; l < r; l++) {
                s120 += psi(k, l) * xi(l, r, L);
                // note: s121 is just s120 w/ derivative wrt (ij) applied to 
                // the first term
                s121 += dpsi(k, l, i, j, psi, params) * xi(l, r, L);
            }

            s10 = psi(k,s) * dpsi(k, r, i, j, psi, params) + 
                psi(k,r) * dpsi(k, s, i, j, psi, params);
            s11 = dpsi(k, r, i, j, psi, params) * s110 + psi(k,r) * s111;
            s12 = dpsi(k, s, i, j, psi, params) * s120 + psi(k,s) * s121;
            s13 = s121 * s110 + s120 * s111;

            s1(k) = s10 + s11 + s12 + s13;
        } // end of for loop computing EACH TERM of s1 summation from 0:(r-1)
        out = -s0 - 1 / psi(r,r) * arma::sum(s1);
    } // end if-else()
    return out; 
}

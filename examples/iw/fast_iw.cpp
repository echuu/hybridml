#include <RcppArmadillo.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
float create_loTriag(arma::mat& L, Rcpp::NumericVector& u);
inline float logmultigamma (unsigned int p, float a);
float cov_loglik (Rcpp::NumericVector& u, Rcpp::List& params); 
float cov_logprior(Rcpp::NumericVector& u, Rcpp::List& params); 


// [[Rcpp::export]]
inline float psi(Rcpp::NumericVector& u, Rcpp::List& params) {
	return -cov_loglik(u, params) - cov_logprior(u, params); 
}

float cov_loglik (Rcpp::NumericVector& u, Rcpp::List& params) {
	int N = params["N"];
	int D = params["D"];
	arma::mat S = params["S"];
	
	arma::mat L(D, D, arma::fill::zeros);

	float logDiagL = create_loTriag(L, u);
	arma::mat inv_L = inv(L);
	arma::mat t_inv_L = inv_L.t();

	float loglik = -0.5 * N * D * std::log(2*M_PI) - N * logDiagL - 
					0.5 * arma::trace(t_inv_L * inv_L * S);

	return loglik;
}


float cov_logprior(Rcpp::NumericVector& u, Rcpp::List& params) {
	arma::mat Omega = params["Omega"];
	unsigned int nu = params["nu"];
	unsigned int D = params["D"];
	arma::mat L(D, D, arma::fill::zeros);

	float logDiagL = create_loTriag(L, u);

	
	float logC = 0.5 * nu * log_det(Omega).real() - 0.5 * (nu * D) * log(2) - 
				logmultigamma(D, nu / 2);

	/**** Computing log of Jacobian term ****/
	float jac = 0;
	for (unsigned int i = 0; i < D; i++) { jac += (D + 1 - (i+1)) * log(L(i, i)); }
	float logJacTerm = D * log(2) + jac;
	/**** End Jacobian computation ****/

	arma::mat inv_L = inv(L);
	arma::mat t_inv_L = inv_L.t();
	float logprior = logC - (nu + D + 1) * logDiagL - 0.5 * arma::trace(t_inv_L * inv_L * Omega) + logJacTerm; 

	return logprior;
}

/* Populars lower diagonal of N x N matrix L with values from u,
	params:
		L: N x N arma matrix
		u: numeric vector
	returns:
		sum of log diag values
*/
float create_loTriag(arma::mat& L, Rcpp::NumericVector& u) {
	if (L.n_rows != L.n_cols)
		throw std::invalid_argument("Unequal dimensions");

	unsigned int D = L.n_rows;
	unsigned int k = 0;	
	float logDiagL = 1;

	for (unsigned int j = 0; j < D; j++) {
		logDiagL *= u[k];
		for (unsigned i = j; i < D; i++) {
			L(i, j) = u[k++];	
		}
	}
	return std::log(logDiagL);
}

/* log multivariate gamma function
	params:
*/
inline float logmultigamma (unsigned int p, float a) {
	float f = 0.25 * p * (p - 1) * std::log(M_PI);
	for (unsigned int i = 1; i <= p; i++)
		f += lgamma(a+0.5-(0.5 * i));
	return f;
}

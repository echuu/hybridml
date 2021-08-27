#define RCPP_ARMADILLO_RETURN_COLVEC_AS_VECTOR
#include <RcppArmadillo.h>
#include <cmath>


typedef unsigned int u_int;


// [[Rcpp::depends(RcppArmadillo)]]
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS

arma::vec chol2vec(arma::mat& M, u_int D);
arma::vec grad_cpp(Rcpp::NumericVector& u, arma::mat& psi_mat, Rcpp::List& params);
float dpsi_cpp(u_int i, u_int j, arma::mat& psi_mat, Rcpp::List& params);
float dpsi_rsij(u_int r, u_int s, u_int i, u_int j, arma::mat& psi_mat, arma::mat& G);

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
arma::vec grad_cpp(Rcpp::NumericVector& u, arma::mat& psi_mat,
	Rcpp::List& params) {


	arma::mat G       = params["G"]; // graph G represented as adjacency matrix
	u_int p           = params["p"]; // dimension of the graph G
	arma::vec edgeInd = params["edgeInd"];
	u_int D           = params["D"]; // dimension of parameter space
	// TODO: implement create_psi_mat() function later; for now, we pass it in
	// arma::mat psi_mat = vec2chol(u, p)

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
float dpsi_cpp(u_int i, u_int j, arma::mat& psi_mat, Rcpp::List& params) {

	arma::mat G     = params["G"];    // graph G represented as adjacency matrix
	u_int p         = params["p"];    // dimension of the graph G
	u_int b         = params["b"];    // degrees of freedom
	arma::vec nu_i  = params["nu_i"]; //

	float d_ij; // derivative of psi wrt psi_ij (summation over 
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
					// Rcpp::Rcout << "G[" << r+1 << ", " << s+1 << "] = " << G(r,s) << std::endl;
				}
			} // end loop over s
		} // end loop over r

	} // end if-else

	return d_ij + psi_mat(i,j);

} // end dpsi() function


/* 
dpsi_rsij() : compute derivative d psi_rs / d psi_ij
*/
// [[Rcpp::export]]
float dpsi_rsij(u_int r, u_int s, u_int i, u_int j, 
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
  	if ((i == r) && (j == s) && G(r, s) == 0) { return 0; } // d wrt to non-free: 0

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
  					tmp_sum(k) = -1/psi_mat(r,r) * dpsi_rsij(k, s, i, j, psi_mat, G) * psi_mat(k,r);
  				} else if (psi_mat(k,r) == 0) {
  					tmp_sum(k) = -1/psi_mat(r,r) * dpsi_rsij(k, r, i, j, psi_mat, G) * psi_mat(k,s);
  				} else {

  					if ((i == j) && (r == i) && (G(r,s) == 0)) {
  						tmp_sum(k) = 1/std::pow(psi_mat(r,r), 2) * psi_mat(k,r) * psi_mat(k,s) -
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
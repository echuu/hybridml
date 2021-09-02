#define RCPP_ARMADILLO_RETURN_COLVEC_AS_VECTOR
#include <RcppArmadillo.h>
#include <cmath>


typedef unsigned int u_int;


// [[Rcpp::depends(RcppArmadillo)]]
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS

/** utility functions **/
arma::mat create_psi_mat_cpp(Rcpp::NumericVector& u, Rcpp::List& params);
arma::vec chol2vec(arma::mat& M, u_int D);


float psi_cpp_mat(arma::mat& psi_mat, Rcpp::List& params);

/** gradient related functions **/
arma::vec grad_cpp(arma::vec& u, Rcpp::List& params);
float dpsi_cpp(u_int i, u_int j, arma::mat& psi_mat, Rcpp::List& params);
float dpsi_rsij(u_int r, u_int s, u_int i, u_int j, 
	arma::mat& psi_mat, arma::mat& G);

/** hessian related functions **/
arma::mat hess_cpp(arma::vec& u, Rcpp::List& params);
float d2(u_int i, u_int j, u_int k, u_int l, 
	arma::mat& psi_mat, Rcpp::List& params);
float d2_rs(u_int r, u_int s, u_int i, u_int j, u_int k, u_int l,
	arma::mat& psi_mat, arma::mat& G);


/** ------------------------------------------------------------------------ **/

// [[Rcpp::export]]
arma::vec approx_integral_0(u_int K, arma::mat& psi_df, Rcpp::List& params) {

	u_int D           = params["D"];    // dimension of parameter space

	arma::vec log_terms(K, arma::fill::zeros);

	arma::mat H_k(D, D, arma::fill::zeros);
	arma::mat H_k_inv(D, D, arma::fill::zeros);
	arma::vec lambda_k(D, arma::fill::zeros);
	arma::vec b_k(D, arma::fill::zeros);
	arma::vec m_k(D, arma::fill::zeros);

	// arma::vec u_k(D, psi_df.submat(0, 0, 0, D-1));

	arma::vec u_k = arma::conv_to< arma::vec >::from(psi_df.row(0));

	for (u_int k = 0; k < K; k++) {

		double val = 0;
		double sign;
		//log_det(val, sign, H_k); 


		log_terms(k) = D / 2 * std::log(2 * M_PI) - 0.5 * val;

	} // end for() over k


	return u_k;
} // end approx_integral() function


// [[Rcpp::export]]
arma::vec approx_integral(u_int K, arma::mat& psi_df, arma::mat& bounds,
	Rcpp::List& params) {

	u_int D           = params["D"];    // dimension of parameter space

	arma::vec log_terms(K, arma::fill::zeros);
	arma::vec G_k(K, arma::fill::zeros);

	arma::mat H_k(D, D, arma::fill::zeros);
	arma::mat H_k_inv(D, D, arma::fill::zeros);
	arma::vec lambda_k(D, arma::fill::zeros);
	arma::vec b_k(D, arma::fill::zeros);
	arma::vec m_k(D, arma::fill::zeros);

	arma::vec u_k(D, arma::fill::zeros);

	for (u_int k = 0; k < K; k++) {

		// Rcpp::Rcout<< k << std::endl;

		u_k = arma::conv_to< arma::vec >::from(psi_df.submat(k, 0, k, D-1));
		// Rcpp::Rcout<< u_k << std::endl;
		H_k = hess_cpp(u_k, params);
		H_k_inv = inv(H_k);
		lambda_k = grad_cpp(u_k, params);
		b_k = H_k * u_k - lambda_k;
		m_k = H_k_inv * b_k;


		// TODO: extract the lower and upper bounds of the k-th partition


		double val = 0;
		double sign;
		log_det(val, sign, H_k); 
		// Rcpp::Rcout << val << std::endl;

		// TODO: load the epmgp code into the same directory so that we can use
		// the EP code directly without having to go back into R env
		
		log_terms(k) = D / 2 * std::log(2 * M_PI) - 0.5 * val + 
			arma::dot(lambda_k, u_k) - 
			(0.5 * u_k.t() * H_k * u_k).eval()(0,0) + 
			(0.5 + m_k.t() * H_k * m_k).eval()(0,0);
		
		// float x =  (m_k.t() * H_k * m_k).eval()(0,0);
		// Rcpp::Rcout << x << std::endl;

	} // end for() over k

	// TODO: find log-sum-exp function in arma
	return log_terms;
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
	float x0, tmp1, tmp2;
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
float psi_cpp(arma::vec& u, Rcpp::List& params) {

	u_int p           = params["p"];    // dimension of the graph G
	u_int D           = params["D"];    // dimension of parameter space
	u_int b           = params["b"];    // degrees of freedom
	arma::vec nu_i    = params["nu_i"]; // see p. 329 of Atay (step 2)
	arma::vec b_i     = params["b_i"];  // see p. 329 of Atay (step 2)
	arma::mat P       = params["P"];    // upper cholesky factor of V_n
	arma::mat G       = params["G"];    // graph G 

	arma::mat psi_mat = create_psi_mat_cpp(u, params);

	float psi_u = p * std::log(2);
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
float psi_cpp_mat(arma::mat& psi_mat, Rcpp::List& params) {

	u_int p           = params["p"];    // dimension of the graph G
	u_int D           = params["D"];    // dimension of parameter space
	u_int b           = params["b"];    // degrees of freedom
	arma::vec nu_i    = params["nu_i"]; // see p. 329 of Atay (step 2)
	arma::vec b_i     = params["b_i"];  // see p. 329 of Atay (step 2)
	arma::mat P       = params["P"];    // upper cholesky factor of V_n
	arma::mat G       = params["G"];    // graph G 

	float psi_u = p * std::log(2);
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
					// Rcpp::Rcout << "G[" << r+1 << ", " << s+1 << \
					// 		"] = " << G(r,s) << std::endl;
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
	float tmp;

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
float d2(u_int i, u_int j, u_int k, u_int l, 
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
			
			tmp(n) = dpsi_rsij(r, s, k, l, psi_mat, G) * 
				dpsi_rsij(r, s, i, j, psi_mat, G) + 
				psi_mat(r,s) * d2_rs(r, s, i, j, k, l, psi_mat, G);
			// tmp(n) = d2_rs(r, s, i, j, k, l, psi_mat, G);
		} // end of if-else

	} // end for() over nonfree elements
	return arma::sum(tmp);
} // end d2() function


// [[Rcpp::export]]
float d2_rs(u_int r, u_int s, u_int i, u_int j, u_int k, u_int l,
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
	if ((r == i) && (r == k) && ((s < j) || (s < l))) { return 0; } 


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

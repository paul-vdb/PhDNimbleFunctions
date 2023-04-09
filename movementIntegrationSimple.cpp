// Base dependencies
#include <iostream>
#include <math.h>
#include <numeric>
#include <vector>
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
double integrateTimeCPP(const double time,
					const double timePrev,
					const arma::mat& traps, 
					const arma::vec& XPrev,
					const arma::vec& S,
					const double lambda,
					const double sigma,	// limiting sigma of OU.
					const double beta,
					const double delta,
					const double toSkip) {
	
	int J = traps.n_rows;

	double tDiff = time - timePrev;
	double N = ceil( tDiff/delta ) + 1.0;	// Make sure it is always at least two points.
	double deltaNew = tDiff / N;
	double ti = deltaNew / 2.0;
	
	double Hk = 0.0;
	double expBT = 0.0;
	double sigmasq2 = pow(sigma, 2)*2;
	double d2 = 0.0;
	
	double x1 = 0.0;
	double x2 = 0.0;
	double new2Var = 0.0;
	
	arma::vec hjLimit = arma::zeros(J);
	for (int j = 0; j < J; j++) {
		hjLimit(j) = exp(-(pow(traps(j,0) - S(0), 2) + pow(traps(j,1) - S(1), 2))/sigmasq2);
	}	
	double hLimit = sum(hjLimit);
	
	for (int i = 0; i < N; i++) {
		if(ti < toSkip/beta){
			expBT = exp(-beta*ti);
			new2Var = (1-exp(-2*beta*ti))*sigmasq2;
			x1 = S(0) + (XPrev(0) - S(0))*expBT;
			x2 = S(1) + (XPrev(1) - S(1))*expBT;
			for (int j = 0; j < J; j++) {
				if(hjLimit(j) > 0.0000001){
					d2 = pow(traps(j,0) - x1, 2) + pow(traps(j,1) - x2, 2);
					Hk += exp(-d2/new2Var);
				}
			}
			ti += deltaNew;
		}else {
			Hk += hLimit;
		}
	}
	Hk *= lambda*deltaNew;
	return Hk;
}
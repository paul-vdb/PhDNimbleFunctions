// Base dependencies
#include <iostream>
#include <math.h>
#include <numeric>
#include <vector>
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
double calculate_hk(const double& t,
					const arma::mat& traps, 
					const arma::vec& X0,
					const arma::vec& S,
					const double& lambda,
					const double& sigma,
					const double& beta)
{	
	int J = traps.n_rows;
	double Hk = 0;
	double sigma2 = pow(sigma, 2);

	double denom = ((1-exp(-2*beta*t))*sigma2);
	for(int j = 0; j < J; j++){
		double d2j = pow(traps(j, 0) - (S(0) + (X0(0) - S(0))*exp(-beta*t)), 2) + pow(traps(j, 1) - (S(1) + (X0(1) - S(1))*exp(-beta*t)), 2);
		Hk += lambda*exp(-d2j/denom);
	}
	return Hk;
}


// [[Rcpp::export]]
double calculateH(const double& time,
					const double& timePrev,
					const arma::mat& traps, 
					const arma::vec& XPrev,
					const arma::vec& S,
					const double& lambda,
					const double& sigma,
					const double& beta,
					const double& delta)
{	
	double tDiff = time - timePrev;
	int J = traps.n_rows;
	//double N = ceil( tDiff/delta );
	//double deltaNew = tDiff / N;
	//double start = delta / 2.0;
	double Hk = 0;
	double sigma2 = pow(sigma, 2);

	//for (int i = 0; i < N; i++) {
		//double ti = start + deltaNew * i;
		double denom = ((1-exp(-2*beta*tDiff))*sigma2);
		for(int j = 0; j < J; j++){
			double d2j = pow(traps(j, 0) - (S(0) + (XPrev(0) - S(0))*exp(-beta*tDiff)), 2) + pow(traps(j, 1) - (S(1) + (XPrev(1) - S(1))*exp(-beta*tDiff)), 2);
			Rcout << "The value of d2j : " << d2j << "\n";			
			//Hk += sum(lambda*exp(-d2j/((1-exp(-2*beta*ti))*sigma2/beta)))*deltaNew;
			Hk += lambda*exp(-d2j/denom);
		}
	//}	
	return Hk;
}

// [[Rcpp::export]]
double integrateTimeCPP_old(const double& time,
					const double& timePrev,
					const arma::mat& traps, 
					const arma::vec& XPrev,
					const arma::vec& S,
					const double& lambda,
					const double& sigma,	// limiting sigma of OU.
					const double& beta,
					const double& delta)
{	
	double tDiff = time - timePrev;
	double N = ceil( tDiff/delta );
	double deltaNew = tDiff / N;
	double start = deltaNew / 2.0;
	double Hk = 0;
	double sigma2 = pow(sigma, 2);

	for (int i = 0; i < N; i++) {
		double ti = start + deltaNew * i;
		arma::vec d2j = pow(traps.col(0) - (S(0) + (XPrev(0) - S(0))*exp(-beta*ti)),2) + pow(traps.col(1)- (S(1) + (XPrev(1) - S(1))*exp(-beta*ti)),2);
		Hk += sum(lambda*exp(-d2j/(2*(1-exp(-2*beta*ti))*sigma2)))*deltaNew;
	}	
	return Hk;
}

// [[Rcpp::export]]
double integrateTimeCPPArma(const double& time,
					const double& timePrev,
					const arma::mat& traps, 
					const arma::vec& XPrev,
					const arma::vec& S,
					const double& lambda,
					const double& sigma,	// limiting sigma of OU.
					const double& beta,
					const double& delta)
{	
	double tDiff = time - timePrev;
	double N = ceil( tDiff/delta );
	double deltaNew = tDiff / N;
	double start = deltaNew / 2.0;
	double Hk = 0;
	double sigma2 = pow(sigma, 2);
	arma::vec traps_sum2 = pow(traps.col(0), 2) + pow(traps.col(1), 2);
	double x1 = 0;
	double x2 = 0;
	double newVar = 0;
	
	for (int i = 0; i < N; i++) {
		double ti = start + deltaNew * i;
		x1 = S(0) + (XPrev(0) - S(0))*exp(-beta*ti);
		x2 = S(1) + (XPrev(1) - S(1))*exp(-beta*ti);
		newVar = (1-exp(-2*beta*ti))*sigma2;
		arma::vec d2j = traps_sum2 - 2*traps.col(0)*x1 - 2*traps.col(1)*x2 + pow(x1, 2) + pow(x2, 2);
		Hk += sum(exp(-d2j/(2*newVar)));
	}
	Hk *= lambda*deltaNew;
	return Hk;
}

// [[Rcpp::export]]
double integrateTimeCPP_basic(const double& time,
					const double& timePrev,
					const arma::mat& traps, 
					const arma::vec& XPrev,
					const arma::vec& S,
					const double& lambda,
					const double& sigma,	// limiting sigma of OU.
					const double& beta,
					const double& delta) {
	
	int J = traps.n_rows;
	double tDiff = time - timePrev;
	double N = ceil( tDiff/delta );
	double deltaNew = tDiff / N;
	double ti = deltaNew / 2.0;
	
	double x1 = 0.0;
	double x2 = 0.0;
	double new2Var = 0.0;
	
	double Hk = 0.0;
	
	for (int i = 0; i < N; i++) {
		x1 = S(0) + (XPrev(0) - S(0))*exp(-beta*ti);
		x2 = S(1) + (XPrev(1) - S(1))*exp(-beta*ti);
		new2Var = 2*(1-exp(-2*beta*ti))*(sigma*sigma);
		for (int j = 0; j < J; j++) {
			Hk += exp(-(pow(traps(j,0) - x1, 2) + pow(traps(j,1) - x2, 2))/new2Var);
		}
		ti += deltaNew;
	}
	Hk *= lambda*deltaNew;
	return Hk;
}

// [[Rcpp::export]]
double integrateTimeCPP(const double& time,
					const double& timePrev,
					const arma::mat& traps, 
					const arma::vec& XPrev,
					const arma::vec& S,
					const double& lambda,
					const double& sigma,	// limiting sigma of OU.
					const double& beta,
					const double& delta,
					const double& toSkip) {
	
	int J = traps.n_rows;
	arma::vec skip = arma::zeros(J);
	if(toSkip == 1){
		for(int k=0; k < J; k++)
		{
			double tmp = (pow(traps(k,0) - S(0), 2) + pow(traps(k,1) - S(1), 2));
			if(tmp > 16*sigma*sigma) { // 4 sigma buffer on distance.
				skip(k) = 1;
			}
		}
	}
	double tDiff = time - timePrev;
	double N = ceil( tDiff/delta );
	double deltaNew = tDiff / N;
	double ti = deltaNew / 2.0;
	
	double x1 = 0.0;
	double x2 = 0.0;
	double new2Var = 0.0;
	
	double Hk = 0.0;
	double stepsDelta = 0.0;
	
	for (int i = 0; i < N; i++) {
		if(ti < 4/beta){
			x1 = S(0) + (XPrev(0) - S(0))*exp(-beta*ti);
			x2 = S(1) + (XPrev(1) - S(1))*exp(-beta*ti);
			new2Var = 2*(1-exp(-2*beta*ti))*(sigma*sigma);
			for (int j = 0; j < J; j++) {
				if(skip(j) == 0){
						Hk += exp(-(pow(traps(j,0) - x1, 2) + pow(traps(j,1) - x2, 2))/new2Var);
				}
			}
			ti += deltaNew;
		}else{
			stepsDelta ++;
		}
	}
	if(stepsDelta > 0)
	{
		for (int j = 0; j < J; j++) {
			if(skip(j) == 0){
					Hk += stepsDelta*exp(-(pow(traps(j,0) - S(0), 2) + pow(traps(j,1) - S(1), 2))/(2*sigma*sigma));
			}
		}
	}
	Hk *= lambda*deltaNew;
	return Hk;
}


double integrateTimeCPP2(const double time,
					const double timePrev,
					const arma::mat& traps, 
					const arma::vec& XPrev,
					const arma::vec& S,
					const double lambda,
					const double sigma,	// limiting sigma of OU.
					const double beta,
					const double delta) {
	
	int J = traps.n_rows;

	double tDiff = time - timePrev;
	double N = ceil( tDiff/delta );
	double deltaNew = tDiff / N;
	double ti = deltaNew / 2.0;
	
	double Hk = 0.0;
	
	double x1 = 0.0;
	double x2 = 0.0;
	double new2Var = 0.0;
	double corr = 0;
	double corr2 = 0;
	double sigma2 = (sigma*sigma)*2;
	
	arma::vec d2Adjusted= arma::zeros(J);
	double stepsDelta = 0.0;
	
	for (int i = 0; i < N; i++) {
		corr = exp(-beta*ti);
		corr2 = pow(corr,2);
		x1 = S(0) + (XPrev(0) - S(0))*corr;
		x2 = S(1) + (XPrev(1) - S(1))*corr;
		new2Var = (1-corr2)*sigma2;
		d2Adjusted = pow(traps.col(0) - x1, 2)/new2Var + pow(traps.col(1) - x2, 2)/new2Var;
		Hk += sum(exp(-d2Adjusted));
		ti += deltaNew;
	}
	
	Hk *= lambda*deltaNew;
	return Hk;
}


// [[Rcpp::export]]
double integrateTimeHackCPP(const double& time,
					const double& timePrev,
					const arma::mat& traps, 
					const arma::vec& XPrev,
					const arma::vec& S,
					const double& lambda,
					const double& sigma,	// limiting sigma of OU.
					const double& beta,
					const double& delta,
					const double& toSkip) {
	
	int J = traps.n_rows;
	arma::vec skip = arma::zeros(J);
	arma::vec Hk0 = arma::zeros(J);	
	arma::vec wgts = arma::zeros(J);

	for(int j=0; j < J; j++)
	{
		Hk0[j] += lambda*exp(-(pow(traps(j,0) - S(0), 2) + pow(traps(j,1) - S(1), 2))/(2*sigma*sigma));
		if(Hk0[j] < 0.000001) { // 4 sigma buffer on distance.
			skip(j) = 1;
		}
	}

	double tDiff = time - timePrev;
	double N = ceil( tDiff/delta );
	double deltaNew = tDiff / N;
	double ti = deltaNew / 2.0;
	
	double x1 = 0.0;
	double x2 = 0.0;
	double new2Var = 0.0;

	double Hk = 0.0;
	double HkSum = sum(Hk0);
	double stepsDelta = 0.0;

	for (int i = 0; i < N; i++) {
		if(ti < 4/beta){
			x1 = S(0) + (XPrev(0) - S(0))*exp(-beta*ti);
			x2 = S(1) + (XPrev(1) - S(1))*exp(-beta*ti);
			new2Var = 2*(1-exp(-2*beta*ti))*(sigma*sigma);
			double wgtsSum = 0;
			for (int j = 0; j < J; j++) {
				if(skip(j) == 0){
					wgts(j) = exp(-(pow(traps(j,0) - x1, 2) + pow(traps(j,1) - x2, 2))/new2Var);
					wgtsSum +=	wgts(j);
				}
				Hk += sum(Hk0%wgts)/wgtsSum;
			}
			ti += deltaNew;
		}else{
			stepsDelta ++;
		}
	}
	if(stepsDelta > 0)
	{
		Hk += arma::dot(Hk0,Hk0)/HkSum;
	}
	Hk *= deltaNew;
	return Hk;
}


// [[Rcpp::export]]
arma::vec integrateTimeCPPMask(const double& time,
					const double& timePrev,
					const arma::mat& traps, 
					const arma::vec& XPrev,
					const arma::vec& S,
					const double& lambda,
					const double& sigma,	// limiting sigma of OU.
					const double& beta,
					const double& delta) {
	
	int J = traps.n_rows;
	double tDiff = time - timePrev;
	double N = ceil( tDiff/delta );
	double deltaNew = tDiff / N;
	double ti = deltaNew / 2.0;
	
	double x1 = 0.0;
	double x2 = 0.0;
	double new2Var = 0.0;
	
	arma::vec Hk = arma::zeros(J);
	
	for (int i = 0; i < N; i++) {
		x1 = S(0) + (XPrev(0) - S(0))*exp(-beta*ti);
		x2 = S(1) + (XPrev(1) - S(1))*exp(-beta*ti);
		new2Var = 2*(1-exp(-2*beta*ti))*(sigma*sigma);
		for (int j = 0; j < J; j++) {
			Hk(j) += exp(-(pow(traps(j,0) - x1, 2) + pow(traps(j,1) - x2, 2))/new2Var);
		}
		ti += deltaNew;
	}
	Hk *= lambda*deltaNew;
	return Hk;
}

// [[Rcpp::export]]
double ESAOUCPP(const double& time,
					const double& timePrev,
					const arma::mat& traps,
					const arma::mat& maskS,
					const arma::mat& maskX,
					const arma::vec& probX,
					const double& areaS,
					const double& areaX,
					const double& lambda,
					const double& sigma,	// limiting sigma of OU.
					const double& beta,
					const double& delta) {
	
	int J = traps.n_rows;
	int Ns = maskS.n_rows;
	int Nx = maskX.n_rows;
	
	double tDiff = time - timePrev;
	double N = ceil( tDiff/delta );
	double deltaNew = tDiff / N;
	double ti = deltaNew / 2.0;

	double ESA = 0;
	
	double x1 = 0.0;
	double x2 = 0.0;
	double new2Var = 0.0;

	for(int si = 0; si < Ns; si++){
		for(int xi = 0; xi < Nx; xi++){
			double Hk = 0;
			for (int i = 0; i < N; i++) {
				x1 = maskS(si,0) + maskX(xi, 0)*exp(-beta*ti);
				x2 = maskS(si,1) + maskX(xi, 1)*exp(-beta*ti);
				new2Var = 2*(1-exp(-2*beta*ti))*(sigma*sigma);
				for (int j = 0; j < J; j++) {
					Hk += exp(-(pow(traps(j,0) - x1, 2) + pow(traps(j,1) - x2, 2))/new2Var);
				}
				ti += deltaNew;
			}
			Hk *= lambda*deltaNew;
			ESA += (1-exp(-Hk))*probX(xi);
		}
	}
	ESA *= areaS*areaX;
	return ESA;
}

// [[Rcpp::export]]
double ESAOUCPPs_x(const double& time,
					const double& timePrev,
					const arma::mat& traps,
					const arma::mat& maskS,
					const double& areaS,
					const double& lambda,
					const double& sigma,	// limiting sigma of OU.
					const double& beta,
					const double& delta) {
	
	int J = traps.n_rows;
	int Ns = maskS.n_rows;
	
	double tDiff = time - timePrev;
	double N = ceil( tDiff/delta );
	double deltaNew = tDiff / N;
	double ti = deltaNew / 2.0;

	double ESA = 0;
	
	double new2Var = 0.0;

	for(int si = 0; si < Ns; si++){
		double Hk = 0;
		for (int i = 0; i < N; i++) {
			new2Var = 2*(1-exp(-2*beta*ti))*(sigma*sigma);
			for (int j = 0; j < J; j++) {
				Hk += exp(-(pow(traps(j,0) - maskS(si,0), 2) + pow(traps(j,1) - maskS(si,1), 2))/new2Var);
			}
			ti += deltaNew;
		}
		Hk *= lambda*deltaNew;
		ESA += (1-exp(-Hk));
	}
	ESA *= areaS;
	return ESA;
}

// [[Rcpp::export]]
double integrateTimeCPP_ARCHIVE(const double& time,
					const double& timePrev,
					const arma::mat& traps, 
					const arma::vec& XPrev,
					const arma::vec& S,
					const double& lambda,
					const double& sigma,	// limiting sigma of OU.
					const double& beta,
					const double& delta,
					const double& toSkip) {
	
	int J = traps.n_rows;
	arma::vec skip = arma::zeros(J);
	if(toSkip == 1){
		for(int k; k < J; k++)
		{
			double tmp = (pow(traps(k,0) - S(0), 2) + pow(traps(k,1) - S(1), 2));
			if(tmp > 16*sigma*sigma) { // 4 sigma buffer on distance.
				skip(k) = 1;
			} 
		}
	}
	double tDiff = time - timePrev;
	double N = ceil( tDiff/delta );
	double deltaNew = tDiff / N;
	double ti = deltaNew / 2.0;
	
	double x1 = 0.0;
	double x2 = 0.0;
	double new2Var = 0.0;
	
	double Hk = 0.0;
	
	for (int i = 0; i < N; i++) {
		x1 = S(0) + (XPrev(0) - S(0))*exp(-beta*ti);
		x2 = S(1) + (XPrev(1) - S(1))*exp(-beta*ti);
		new2Var = 2*(1-exp(-2*beta*ti))*(sigma*sigma);
		for (int j = 0; j < J; j++) {
			if(skip(j) == 0){
					Hk += exp(-(pow(traps(j,0) - x1, 2) + pow(traps(j,1) - x2, 2))/new2Var);
			}
		}
		ti += deltaNew;
	}
	Hk *= lambda*deltaNew;
	return Hk;
}

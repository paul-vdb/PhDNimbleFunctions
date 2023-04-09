// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <RcppNumerical.h>
using namespace Numer;
using Eigen::Map; 

class OUHazardAll: public Func
{
private:
    const MapMat X;
    const MapVec S;	
    const MapVec X0;
	const double sigma;
	const double beta;
	const double lambda;
public:
    OUHazardAll(const MapMat x_, const MapMat s_, const MapVec x0_,
		const double& sigma_, const double& beta_, const double& lambda_): 
		X(x_), S(s_), X0(x0_),
		sigma(sigma_), beta(beta_), lambda(lambda_) {}

    double operator()(const double& t) const
    {
		Eigen::MatrixXf D = X.rowwise() -  (S + (X0-S)*exp(-beta*t));
		double sigmat2 = 2*(1-exp(-2*beta*t))*pow(sigma,2);
		Eigen::VectorXd eD = D.square().rowwise().sum()
        return lambda*exp(-(pow(d1,2) + pow(d2,2))/sigmat2);
    }
};

// [[Rcpp::export]]
Rcpp::List integrate_test(const double x1, const double x2, 
	const double s1, const double s2,
	const double x01, const double x02, 
	double sigma, double lambda, double beta,
	double lower, double upper)
{

    OUHazard f(x1, x2, s1, s2, x01, x02, sigma, beta, lambda);
    double err_est;
    int err_code;
    const double res = integrate(f, lower, upper, err_est, err_code);
    return Rcpp::List::create(
        Rcpp::Named("approximate") = res,
        Rcpp::Named("error_estimate") = err_est,
        Rcpp::Named("error_code") = err_code
    );
}

// [[Rcpp::export]]
double integrate_OU(const double& x1, const double& x2, 
	const double& s1, const double& s2,
	const double& x01, const double& x02, 
	const double& sigma, double& lambda, double& beta,
	double& lower, double& upper )
{

    OUHazard f(x1, x2, s1, s2, x01, x02, sigma, beta, lambda);
    double err_est;
    int err_code;
    const double res = integrate(f, lower, upper, err_est, err_code);
    return res;
}

//// [[Rcpp::export]]
// double cumulativeIntensityCpp(const double& time,
					// const double& timePrev,
					// const arma::mat& traps, 
					// const arma::vec& XPrev,
					// const arma::vec& S,
					// const double& lambda,
					// const double& sigma,
					// const double& beta)
// {
	// double H = 0;
	// int J = traps.n_rows;
	// double tDiff = time - timePrev;
	// double s1 = S(0);
	// double s2 = S(1);
	// double x01 = XPrev(0);
	// double x02 = XPrev(1);
	// double low = 0;
	
	//for (int j = 0; j < J; j++){
		// double x1 = traps(j,0);
		// double x2 = traps(j,1);
		// H += integrate_OU(x1 = x1, x2 = x2, 
			// s1 = s1, s2 = s2, 
			// x01 = x01, x02 = x02, 
			// sigma = sigma, beta = beta, lambda = lambda,
			// lower = low, upper = tDiff);
	//}
	// return H;
// }
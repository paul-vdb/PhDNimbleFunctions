#include <RcppArmadillo.h>
#include <RcppNumerical.h>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Numer;

class OUHazard: public Func
{
private:
    double X1;
    double X2;
	double S1;
	double S2;
	double X01;
	double X02;
	double Sigma;
	double Beta;
	double Lambda;
public:
    OUHazard(const double x1_, const double x2_,
		const double s1_, const double s2_,
		const double x01_, const double x02_,
		const double sigma_, const double beta_, const double lambda_): 
		X1(x1_), X2(x2_), 
		S1(s1_), S2(s2_), 
		X01(x01_), X02(x02_), 
		Sigma(sigma_), Beta(beta_), Lambda(lambda_) {}

    double operator()(const double& t) const
    {
		double d1 = (X1 - (S1 + (X01-S1)*exp(-Beta*t)));
		double d2 = (X2 - (S2 + (X02-S2)*exp(-Beta*t)));
		double sigmat2 = 2*(1-exp(-2*Beta*t))*pow(Sigma,2);
        return Lambda*exp(-(pow(d1,2) + pow(d2,2))/sigmat2);
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
double integrate_OU(const double x1, const double x2, 
	const double s1, const double s2,
	const double x01, const double x02, 
	const double sigma, double lambda, double beta,
	double lower, double upper )
{

    OUHazard f(x1, x2, s1, s2, x01, x02, sigma, beta, lambda);
    double err_est;
    int err_code;
    const double res = integrate(f, lower, upper, err_est, err_code);
    return res;
}

//// [[Rcpp::export]]
// double cumulativeIntensityCpp(const double ti,
					// const double timePrev,
					// const arma::mat traps, 
					// const arma::vec XPrev,
					// const arma::vec S,
					// const double lambda,
					// const double sigma,
					// const double beta)
// {
	// double H = 0;
	// int J = traps.n_rows;
	// double tDiff = ti - timePrev;
	// double s1 = S(0);
	// double s2 = S(1);
	// double x01 = XPrev(0);
	// double x02 = XPrev(1);
	// double low = 0;
	// int j = 0; // PAUL YOU IDIOT REMOVE THIS.
	// for (int j = 0; j < J; j++){
		// double x1 = traps(j,0);
		// double x2 = traps(j,1);
		// H += integrate_OU(x1 = x1, x2 = x2, 
			// s1 = s1, s2 = s2, 
			// x01 = x01, x02 = x02, 
			// sigma = sigma, beta = beta, lambda = lambda,
			// lower = low, upper = tDiff);
	// }
	// return H;
// }
#ifndef GAUSSIAN2D_H
#define GAUSSIAN2D_H

#include "RNG.h"
#include "Eigen/Dense"

class Gaussian2D
{
	private:
		double mean_1_, mean_2_, rho_, sigma_1_, sigma_2_, log_Norm;
		// Covariance matrix
		Eigen::Matrix2d Cov;
		// Inverse covariance matrix
		Eigen::Matrix2d Cov_inv;
		// Choleski decomposition of Sigma (upper tri: Mt M = Cov)
		Eigen::Matrix2d M;
		
	public:
		Gaussian2D(double mean_1, double mean_2, double rho, double sigma_1, double sigma_2);
		std::vector<double> generate(DNest4::RNG& rng) const;
		double perturb(std::vector<double>& x, DNest4::RNG& rng) const;
		
//		virtual void setpars(double) {};
//		virtual void setpars(double, double) {};
		
		std::vector<double> cdf(std::vector<double> x) const;
		// It is not CDF actually
		std::vector<double> cdf_inverse(std::vector<double> u) const;
		double log_pdf(std::vector<double> x) const;
};

#endif //GAUSSIAN2D_H

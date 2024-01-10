#include <stdexcept>
#include <Distributions/Gaussian.h>
#include <iostream>
#include "Gaussian2D.h"


Gaussian2D::Gaussian2D(std::vector<double> mean, double rho, double sigma_1, double sigma_2) :
mean_(mean.data()), rho_(rho), sigma_1_(sigma_1), sigma_2_(sigma_2)
{
	if(sigma_1 <=0 || sigma_2 <=0)
		throw std::domain_error("Gaussian2D distribution must have positive widths.");
	if(abs(rho) >= 1)
		throw std::domain_error("Gaussian2D distribution must have correlation > 1.");
		
	Cov << sigma_1_*sigma_1_, rho_*sigma_1_*sigma_2_,
	       rho_*sigma_1_*sigma_2_, sigma_2_*sigma_2_;

	Cov_inv << 1./(sigma_1_*sigma_1_), -rho_/(sigma_1_*sigma_2_),
		       -rho_/(sigma_1_*sigma_2_), 1./(sigma_2_*sigma_2_);
	Cov_inv *= (1./(1. - rho_*rho_));
	
	M << sigma_1_, rho_*sigma_2_,
	     0., sigma_2_*sqrt(1. - rho_*rho_);
			
	log_Norm = -log(2.*M_PI*sqrt(Cov.determinant()));
}


std::vector<double> Gaussian2D::cdf(std::vector<double> x) const
{
	Eigen::Vector2d xx(x.data());
	Eigen::Vector2d diff = xx - mean_;
	Eigen::Vector2d x_st = M.inverse()*diff;
	DNest4::Gaussian gaussian(0, 1);
	return {gaussian.cdf(x_st[0]), gaussian.cdf(x_st[1])};
	
}
std::vector<double> Gaussian2D::cdf_inverse(std::vector<double> u) const
{
	// First, generate two cdf_inverse for standard normal
	DNest4::Gaussian gaussian(0, 1);
	Eigen::Vector2d u_st = {gaussian.cdf_inverse(u[0]), gaussian.cdf_inverse(u[1])};
	Eigen::Vector2d sample = u_st.transpose()*M;
	sample += mean_;
	return {sample.data(), sample.data() + sample.size()};
}

double Gaussian2D::log_pdf(std::vector<double> x) const
{
	Eigen::Vector2d xx(x.data());
	Eigen::Vector2d diff = xx - mean_;
	
	return log_Norm - 0.5*diff.transpose()*Cov_inv*diff;
}

std::vector<double> Gaussian2D::generate(DNest4::RNG& rng) const
{
	Eigen::Vector2d x_st(rng.randn(), rng.randn());
    Eigen::Vector2d sample = x_st.transpose()*M;
	sample += mean_;
	return {sample.data(), sample.data() + sample.size()};
}

double Gaussian2D::perturb(std::vector<double>& x, DNest4::RNG& rng) const
{
	std::vector<double> xx;
	xx.resize(2);
    xx = cdf(x);
	// This works
//	std::cout << "Before adding randh : " << xx[0] << "\n";
	xx[0] += rng.randh();
	xx[1] += rng.randh();
//	std::cout << "After adding randh: " << xx[0] << "\n";
	// ?
//	std::cout << "Before : " << xx[0] << "\n";
	DNest4::wrap(xx[0], 0.0, 1.0);
	DNest4::wrap(xx[1], 0.0, 1.0);
//	std::cout << "After : " << xx[0] << "\n";
    xx = cdf_inverse(xx);
	x = xx;
	
    return 0.0;
}

#include <stdexcept>
#include <Distributions/Gaussian.h>
#include <iostream>
#include "Gaussian2D.h"


Gaussian2D::Gaussian2D(double mean_1, double mean_2, double rho, double sigma_1, double sigma_2) :
mean_1_(mean_1), mean_2_(mean_2), rho_(rho), sigma_1_(sigma_1), sigma_2_(sigma_2)
{
	log_Norm = -log(2.*M_PI*sqrt(1. - rho*rho)*sigma_1_*sigma_2_);
}


std::vector<double> Gaussian2D::cdf(std::vector<double> x) const
{
	DNest4::Gaussian norm_1(mean_1_, sigma_1_);
	double u_1 = norm_1.cdf(x[0]);
	DNest4::Gaussian norm_2(mean_2_ + rho_*sigma_2_/sigma_1_*(x[0] - mean_1_), sqrt(1. - rho_*rho_)*sigma_2_);
	double u_2 = norm_2.cdf(x[1]);
	return {u_1, u_2};
	
}
std::vector<double> Gaussian2D::cdf_inverse(std::vector<double> u) const
{
	DNest4::Gaussian norm_1(mean_1_, sigma_1_);
	double x_1 = norm_1.cdf_inverse(u[0]);
	DNest4::Gaussian norm_2(mean_2_ + rho_*sigma_2_/sigma_1_*(x_1 - mean_1_), sqrt(1. - rho_*rho_)*sigma_2_);
	double x_2 = norm_2.cdf_inverse(u[1]);
	return {x_1, x_2};
}

double Gaussian2D::log_pdf(std::vector<double> x) const
{
	return log_Norm - 0.5*(x[0]*x[0]/(sigma_1_*sigma_1_) - 2.*rho_*x[0]*x[1]/(sigma_1_*sigma_2_) + x[1]*x[1]/(sigma_2_*sigma_2_))/(1. - rho_*rho_);
}

std::vector<double> Gaussian2D::generate(DNest4::RNG& rng) const
{
	return cdf_inverse({rng.rand(), rng.rand()});
}

double Gaussian2D::perturb(std::vector<double>& x, DNest4::RNG& rng) const
{
	x = cdf(x);
	x[0] += rng.randh();
	DNest4::wrap(x[0], 0.0, 1.0);
	x[1] += rng.randh();
	DNest4::wrap(x[1], 0.0, 1.0);
	x = cdf_inverse(x);
	
	return 0.0;
}

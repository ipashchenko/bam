#include "Component.h"
#include "Data.h"
#include <cmath>
#include <iostream>
#include <RNG.h>


void Component::ft(double nu, const std::valarray<double>& u, const std::valarray<double>& v)
{
	switch (shape)
	{
	case Gaussian:
		{
			std::valarray<double> theta;
			double c;
			std::valarray<double> ft;
			
			auto position = get_pos(nu);
			double logsize = get_logsize(nu);
			double logflux = get_logflux(nu);
			// Phase shift due to not being in a phase center
			theta = 2*M_PI*mas_to_rad*(u*position.first + v*position.second);
			// Calculate FT of a Gaussian in a phase center
			c = pow(M_PI*exp(logsize)*mas_to_rad, 2)/(4.*log(2.));
			// TODO: Keep u*u and v*v already computed
			ft = exp(logflux - c*(u*u + v*v));
			
			// Prediction of visibilities
			mu_real = ft*cos(theta);
			mu_imag = ft*sin(theta);
		}
	case Sphere:
		{
			std::valarray<double> theta;
			std::valarray<double> ft;
			
			auto position = get_pos(nu);
			double logD = get_logsize(nu);
			double logflux = get_logflux(nu);
			// Phase shift due to not being in a phase center
			theta = 2*M_PI*mas_to_rad*(u*position.first + v*position.second);
			// Calculate FT of a Sphere in a phase center
			
			std::valarray<double> pi_D_rho = M_PI*exp(logD)*sqrt(u*u + v*v);
			ft = 3*exp(logflux)*(sin(pi_D_rho) - pi_D_rho*cos(pi_D_rho))/pow(pi_D_rho, 3);
			
			// Prediction of visibilities
			mu_real = ft*cos(theta);
			mu_imag = ft*sin(theta);
		}
	}
}


JetComponent::JetComponent(Shape new_shape) : Component() {
	shape = new_shape;
}

JetComponent::JetComponent(const JetComponent &other)
{
    dx_ = other.dx_;
    dy_ = other.dy_;
    logsize_ = other.logsize_;
	lognu_max_ = other.lognu_max_;
	logS_max_ = other.logS_max_;
	alpha_thick_ = other.alpha_thick_;
	alpha_thin_ = other.alpha_thin_;
    mu_real = other.get_mu_real();
    mu_imag = other.get_mu_imag();
	shape = other.shape;
}


JetComponent* JetComponent::clone()
{
    return new JetComponent(*this);
}

double JetComponent::get_logsize(double nu)
{
	return logsize_;
}

double JetComponent::get_logflux(double nu)
{
	return logS_max_ + alpha_thick_*log(nu/exp(lognu_max_)) - log(1 - exp(-1)) + log(1 - exp(-pow(nu/exp(lognu_max_), alpha_thin_-alpha_thick_)));
}

std::pair<double, double> JetComponent::get_pos(double nu)
{
	return {dx_, dy_};
}


void JetComponent::print(std::ostream &out) const
{
	out << dx_ << "\t" << dy_ << "\t" << logsize_ << "\t" << lognu_max_ << "\t" << logS_max_ << "\t" << alpha_thick_ << "\t" << alpha_thin_ << "\t";

}

std::string JetComponent::description() const
{
    std::string descr;
    descr += "dx\tdy\tlogsize\tlognu_max\tlogS_max\talpha_thick\talpha_thin";
    return descr;
}

void JetComponent::from_prior(DNest4::RNG &rng) {
	 DNest4::Gaussian gaussian_pos(0.0, 5.0);
	 DNest4::Gaussian gaussian_logflux(-2.0, 1.0);
	 DNest4::Gaussian gaussian_logsize(-1.0, 2.0);
	 DNest4::Gaussian gaussian_alpha_thick(1.5, 0.5);
	 DNest4::Gaussian gaussian_alpha_thin(-1.0, 0.25);
	 DNest4::Uniform uniform_numax(-3, 3);
     dx_ = gaussian_pos.generate(rng);
     dy_ = gaussian_pos.generate(rng);
     logS_max_ = gaussian_logflux.generate(rng);
     logsize_ = gaussian_logsize.generate(rng);
	 lognu_max_ = uniform_numax.generate(rng);
	 alpha_thick_ = gaussian_alpha_thick.generate(rng);
	 alpha_thin_ = gaussian_alpha_thin.generate(rng);
}

double JetComponent::perturb(DNest4::RNG &rng) {
    double log_H = 0.;
    int which = rng.rand_int(7);
    if(which == 0)
    {
		DNest4::Gaussian gaussian_pos(0.0, 5.0);
		log_H += gaussian_pos.perturb(dx_, rng);
    }
    else if(which == 1)
    {
		DNest4::Gaussian gaussian_pos(0.0, 5.0);
		log_H += gaussian_pos.perturb(dy_, rng);
    }
    else if(which == 2)
    {
		DNest4::Gaussian gaussian_logflux(-2.0, 1.0);
		log_H += gaussian_logflux.perturb(logS_max_, rng);
    }
    else if(which == 3)
    {
		DNest4::Gaussian gaussian_logsize(-1.0, 2.0);
		log_H += gaussian_logsize.perturb(logsize_, rng);
    }
	else if(which == 4)
	{
		DNest4::Uniform uniform_numax(-3, 3);
		log_H += uniform_numax.perturb(lognu_max_, rng);

	}
	else if(which == 5)
	{
		DNest4::Gaussian gaussian_alpha_thick(1.5, 0.5);
		log_H += gaussian_alpha_thick.perturb(alpha_thick_, rng);
	}
	else
	{
		DNest4::Gaussian gaussian_alpha_thin(-1.0, 0.25);
		log_H += gaussian_alpha_thin.perturb(alpha_thin_, rng);
	}
    return log_H;
}

CoreComponent::CoreComponent(Shape new_shape) : Component() {
	shape = new_shape;
}

CoreComponent::CoreComponent(const CoreComponent &other)
{
	a_ = other.a_;
	PA_ = other.PA_;
	logsize_1_ = other.logsize_1_;
	k_r_ = other.k_r_;
//	logS_1_ = other.logS_1_;
//	alpha_ = other.alpha_;
	lognu_max_ = other.lognu_max_;
	logS_max_ = other.logS_max_;
	alpha_thick_ = other.alpha_thick_;
	alpha_thin_ = other.alpha_thin_;
	shape = other.shape;
}

CoreComponent *CoreComponent::clone()
{
	return new CoreComponent(*this);
}

double CoreComponent::get_logsize(double nu)
{
	return logsize_1_ - (1./k_r_)*log(nu);
}

double CoreComponent::get_logflux(double nu)
{
//	return logS_1_ + alpha_*log(nu);
	return logS_max_ + alpha_thick_*log(nu/exp(lognu_max_)) - log(1 - exp(-1)) + log(1 - exp(-pow(nu/exp(lognu_max_), alpha_thin_-alpha_thick_)));
}

std::pair<double, double> CoreComponent::get_pos(double nu)
{
	double distance = a_*pow(nu, -1/k_r_);
	return {distance*sin(PA_), distance*cos(PA_)};
}

void CoreComponent::print(std::ostream &out) const
{
//	out << a_ << "\t" << PA_ << "\t" << logsize_1_ << "\t" << k_r_ << "\t" << logS_1_ << "\t" << alpha_ << "\t";
	out << a_ << "\t" << PA_ << "\t" << logsize_1_ << "\t" << k_r_ << "\t" << lognu_max_ << "\t" << logS_max_ << "\t" << alpha_thick_ << "\t" << alpha_thin_ << "\t";
}

std::string CoreComponent::description() const
{
	std::string descr;
//	descr += "a\tPA\tlogsize_1\tk_r\tlogS_1\talpha";
	descr += "a\tPA\tlogsize_1\tk_r\tlognu_max\tlogS_max\talpha_thick\talpha_thin";
	return descr;
}

void CoreComponent::from_prior(DNest4::RNG &rng)
{
	DNest4::TruncatedCauchy cauchy_pos(0.0, 1.0, 0.0, 10.0);
	DNest4::Uniform gaussian_direction(0.0, 2*M_PI);
	DNest4::Gaussian gaussian_logsize(-1.0, 2.0);
	DNest4::Gaussian gaussian_logflux(-1.0, 1.0);
//	DNest4::Uniform uniform_alpha(-1.5, 1.5);
	DNest4::Uniform uniform_k(0.0, 3.0);
	DNest4::Gaussian gaussian_alpha_thick(1.5, 0.5);
	DNest4::Gaussian gaussian_alpha_thin(-1.0, 0.25);
	DNest4::Uniform uniform_numax(-3, 3);
	a_ = cauchy_pos.generate(rng);
	PA_ = gaussian_direction.generate(rng);
	logsize_1_ = gaussian_logsize.generate(rng);
	k_r_ = uniform_k.generate(rng);
//	logS_1_ = gaussian_logflux.generate(rng);
//	alpha_ = uniform_alpha.generate(rng);
	logS_max_ = gaussian_logflux.generate(rng);
	lognu_max_ = uniform_numax.generate(rng);
	alpha_thick_ = gaussian_alpha_thick.generate(rng);
	alpha_thin_ = gaussian_alpha_thin.generate(rng);
}

double CoreComponent::perturb(DNest4::RNG &rng)
{
	double log_H = 0.;
	int which = rng.rand_int(8);
	if(which == 0)
	{
		DNest4::TruncatedCauchy cauchy_pos(0.0, 1.0, 0.0, 10.0);
		log_H += cauchy_pos.perturb(a_, rng);
	}
	else if(which == 1)
	{
		DNest4::Uniform gaussian_direction(0.0, 2*M_PI);
		log_H += gaussian_direction.perturb(PA_, rng);
	}
	else if(which == 2)
	{
		DNest4::Gaussian gaussian_logsize(-1.0, 2.0);
		log_H += gaussian_logsize.perturb(logsize_1_, rng);
	}
	else if(which == 3)
	{
		DNest4::Gaussian gaussian_logflux(-2.0, 1.0);
		log_H += gaussian_logflux.perturb(logS_max_, rng);
	}
	
//	else if(which == 3)
//	{
//		DNest4::Gaussian gaussian_logflux(-1.0, 1.0);
//		log_H += gaussian_logflux.perturb(logS_1_, rng);
//	}
//	else if(which == 4)
//	{
//		DNest4::Uniform uniform_alpha(-1.5, 1.5);
//		log_H += uniform_alpha.perturb(alpha_, rng);
//	}
	else if(which == 4)
	{
		DNest4::Uniform uniform_k(0.0, 3.0);
		log_H += uniform_k.perturb(k_r_, rng);
	}
	else if(which == 5)
	{
		DNest4::Uniform uniform_numax(-3, 3);
		log_H += uniform_numax.perturb(lognu_max_, rng);
		
	}
	else if(which == 6)
	{
		DNest4::Gaussian gaussian_alpha_thick(1.5, 0.5);
		log_H += gaussian_alpha_thick.perturb(alpha_thick_, rng);
	}
	else
	{
		DNest4::Gaussian gaussian_alpha_thin(-1.0, 0.25);
		log_H += gaussian_alpha_thin.perturb(alpha_thin_, rng);
	}
	return log_H;
}

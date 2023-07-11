#include "Component.h"
#include "Data.h"
#include <cmath>
#include <iostream>
#include <RNG.h>


//EGComponent::EGComponent() : dx_(0.0), dy_(0.0), logbmaj_(0.0), e_(1.0), bpa_(0.0), lognu_max(0.0),
//							 logS_max(1.0), alpha_thick(1.0), alpha_thin(-1.0) {}
//
//
//void EGComponent::ft(double nu, std::valarray<double> u, std::valarray<double> v)
//{
//    std::valarray<double> theta;
//    double c;
//    std::valarray<double> b;
//    std::valarray<double> ft;
//
//    // Phase shift due to not being in a phase center
//    theta = 2*M_PI*mas_to_rad*(u*dx_+v*dy_);
//    // Calculate FT of a Gaussian in a phase center
//    c = pow(M_PI*exp(logbmaj_)*mas_to_rad, 2)/(4.*log(2.));
//    b = e_*e_*pow((u*cos(bpa_) - v*sin(bpa_)), 2) + pow((u*sin(bpa_)+v*cos(bpa_)), 2);
//	double logflux = logS_max + alpha_thick*log(nu/exp(lognu_max)) - log(1 - exp(-1)) + log(1 - exp(-pow(nu/exp(lognu_max), alpha_thin-alpha_thick)));
//    ft = exp(logflux - c*b);
//
//    // Prediction of visibilities
//    mu_real = ft*cos(theta);
//    mu_imag = ft*sin(theta);
//}
//
//
//void EGComponent::print(std::ostream &out) const
//{
//    out << dx_ << "\t" << dy_ << "\t" << logbmaj_ << "\t" << e_ << "\t" << bpa_ << "\t" << lognu_max << "\t" << logS_max << "\t" << alpha_thick << "\t" << alpha_thin << "\n";
//}
//
//
//CGComponent::CGComponent() : EGComponent() {}
//
//
//void CGComponent::ft(double nu, std::valarray<double> u, std::valarray<double> v)
//{
//    std::valarray<double> theta;
//    double c;
//    //std::valarray<double> b;
//    std::valarray<double> ft;
//
//    // Phase shift due to not being in a phase center
//    theta = 2*M_PI*mas_to_rad*(u*dx_+v*dy_);
//    // Calculate FT of a Gaussian in a phase center
//    c = pow(M_PI*exp(logbmaj_)*mas_to_rad, 2)/(4.*log(2.));
//    //b = e_*e_*pow((u*cos(bpa_) - v*sin(bpa_)), 2) + pow((u*sin(bpa_)+v*cos(bpa_)), 2);
//    // TODO: Keep u*u and v*v already computed
//	double logflux = logS_max + alpha_thick*log(nu/exp(lognu_max)) - log(1 - exp(-1)) + log(1 - exp(-pow(nu/exp(lognu_max), alpha_thin-alpha_thick)));
//    ft = exp(logflux - c*(u*u + v*v));
//
//    // Prediction of visibilities
//    mu_real = ft*cos(theta);
//    mu_imag = ft*sin(theta);
//}
//
//
//void CGComponent::print(std::ostream &out) const
//{
//	out << dx_ << "\t" << dy_ << "\t" << logbmaj_ << "\t" << lognu_max << "\t" << logS_max << "\t" << alpha_thick << "\t" << alpha_thin << "\n";
//
//}
//
//
//void CGComponent::from_prior(DNest4::RNG &rng) {
////     const std::valarray<double>& u = Data::get_instance().get_u();
////     std::valarray<double> zero (0.0, u.size());
////     mu_real = zero;
////     mu_imag = zero;
//     // Normal diffuse prior for x & y
//	 // FIXME: If core - use tight prior around phase center (see Kima)
//	 DNest4::Gaussian gaussian_pos(0.0, 5.0);
//	 DNest4::Gaussian gaussian_logflux(-2.0, 1.0);
//	 DNest4::Gaussian gaussian_logsize(-1.0, 2.0);
//	 DNest4::Gaussian gaussian_alpha_thick(1.0, 0.5);
//	 DNest4::Gaussian gaussian_alpha_thin(-1.0, 0.5);
//	 DNest4::Uniform uniform(1, 10);
//     dx_ = gaussian_pos.generate(rng);
//     dy_ = gaussian_pos.generate(rng);
//     // Log-normal prior for flux and bmaj
//     logS_max = gaussian_logflux.generate(rng);
//     logbmaj_ = gaussian_logsize.generate(rng);
//	 lognu_max = uniform.generate(rng);
//	 alpha_thick = gaussian_alpha_thick.generate(rng);
//	 alpha_thin = gaussian_alpha_thin.generate(rng);
//    // Dependent priors - from MOJAVE 15 GHz modelfits.
//    //logflux_ = -1.74 + 1.62*rng.randn();
//    //logbmaj_ = -0.71 - 0.44*logflux_ + 0.84*rng.randn();
//    // Dependent priors - modified for 43 GHz.
////    logflux_ = -2.29 + 1.62*rng.randn();
////    logbmaj_ = -1.71 - 0.44*logflux_ + 0.84*rng.randn();
//}
//
//
//double CGComponent::perturb(DNest4::RNG &rng) {
//    double log_H = 0.;
//    int which = rng.rand_int(7);
//    if(which == 0)
//    {
//		DNest4::Gaussian gaussian_pos(0.0, 5.0);
//		log_H += gaussian_pos.perturb(dx_, rng);
//    }
//    else if(which == 1)
//    {
//		DNest4::Gaussian gaussian_pos(0.0, 5.0);
//		log_H += gaussian_pos.perturb(dy_, rng);
//    }
//    else if(which == 2)
//    {
//		DNest4::Gaussian gaussian_logflux(-2.0, 1.0);
//		log_H += gaussian_logflux.perturb(logS_max, rng);
//
////        log_H -= -0.5*pow((logflux_+1.74)/1.62, 2);
////        logflux_ += 1.62*rng.randh();
////        log_H += -0.5*pow((logflux_+1.74)/1.62, 2);
//    }
//    else if(which == 3)
//    {
//		DNest4::Gaussian gaussian_logsize(-1.0, 2.0);
//		log_H += gaussian_logsize.perturb(logbmaj_, rng);
//
////        log_H -= -0.5*pow((logbmaj_+0.71+0.44*logflux_)/0.84, 2);
////        logbmaj_ += 0.84*rng.randh();
////        log_H += -0.5*pow((logbmaj_+0.71+0.44*logflux_)/0.84, 2);
//    }
//	else if(which == 4)
//	{
//		DNest4::Uniform uniform(1, 10);
//		log_H += uniform.perturb(lognu_max, rng);
//
//	}
//	else if(which == 5)
//	{
//		DNest4::Gaussian gaussian_alpha_thick(1.0, 0.5);
//		log_H += gaussian_alpha_thick.perturb(alpha_thick, rng);
//	}
//	else
//	{
//		DNest4::Gaussian gaussian_alpha_thin(-1.0, 0.5);
//		log_H += gaussian_alpha_thin.perturb(alpha_thin, rng);
//	}
//    return log_H;
//}
//
//
//CGComponent *CGComponent::clone() {
//    return new CGComponent(*this);
//}
//
//
//CGComponent::CGComponent(const CGComponent &other) {
//    dx_ = other.dx_;
//    dy_ = other.dy_;
//    logbmaj_ = other.logbmaj_;
//	lognu_max = other.lognu_max;
//	logS_max = other.logS_max;
//	alpha_thick = other.alpha_thick;
//	alpha_thin = other.alpha_thin;
//    mu_real = other.get_mu_real();
//    mu_imag = other.get_mu_imag();
//}
//
//
//std::string CGComponent::description() const {
//    std::string descr;
//    descr += "x y logbmaj nu_mas logS_max alpha_thick alpha_thin";
//    return descr;
//}
//
//
//CoreCGComponent::CoreCGComponent() {}
//void CoreCGComponent::ft(double nu, std::valarray<double> u, std::valarray<double> v)
//{
//	std::valarray<double> theta;
//	double c;
//	std::valarray<double> ft;
//
//	// Size at given frequency
//	logbmaj_ = log_theta_1 - (1./k_r)*log(nu);
//	// Flux at given frequency
//	double logflux = logS_1 + alpha*log(nu);
//
//	// Phase shift due to not being in a phase center
//	theta = 2*M_PI*mas_to_rad*(u*dx_+v*dy_);
//	// Calculate FT of a Gaussian in a phase center
//	c = pow(M_PI*exp(logbmaj_)*mas_to_rad, 2)/(4.*log(2.));
//	ft = exp(logflux - c*(u*u + v*v));
//
//	// Prediction of visibilities
//	mu_real = ft*cos(theta);
//	mu_imag = ft*sin(theta);
//}
//
//void CoreCGComponent::from_prior(DNest4::RNG &rng)
//{
//	dx_ = 0.0;
//	dy_ = 0.0;
//	logbmaj_ = 0.0;
//}
//
//
//

GaussianComponent::GaussianComponent() {}

void GaussianComponent::ft(double nu, std::valarray<double> u, std::valarray<double> v)
{
	std::valarray<double> theta;
    double c;
    //std::valarray<double> b;
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




JetGaussianComponent::JetGaussianComponent() : GaussianComponent() {}

JetGaussianComponent::JetGaussianComponent(const JetGaussianComponent &other)
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
}

JetGaussianComponent* JetGaussianComponent::clone()
{
    return new JetGaussianComponent(*this);
}

double JetGaussianComponent::get_logsize(double nu)
{
	return logsize_;
}

double JetGaussianComponent::get_logflux(double nu)
{
	return logS_max_ + alpha_thick_*log(nu/exp(lognu_max_)) - log(1 - exp(-1)) + log(1 - exp(-pow(nu/exp(lognu_max_), alpha_thin_-alpha_thick_)));
}

std::pair<double, double> JetGaussianComponent::get_pos(double nu)
{
	return {dx_, dy_};
}


void JetGaussianComponent::print(std::ostream &out) const
{
	out << dx_ << "\t" << dy_ << "\t" << logsize_ << "\t" << lognu_max_ << "\t" << logS_max_ << "\t" << alpha_thick_ << "\t" << alpha_thin_ << "\t";

}

std::string JetGaussianComponent::description() const
{
    std::string descr;
    descr += "dx dy logsize lognu_mas logS_max alpha_thick alpha_thin";
    return descr;
}

void JetGaussianComponent::from_prior(DNest4::RNG &rng) {
	 DNest4::Gaussian gaussian_pos(0.0, 5.0);
	 DNest4::Gaussian gaussian_logflux(-2.0, 1.0);
	 DNest4::Gaussian gaussian_logsize(-1.0, 2.0);
	 DNest4::Gaussian gaussian_alpha_thick(1.0, 0.5);
	 DNest4::Gaussian gaussian_alpha_thin(-1.0, 0.5);
	 DNest4::Uniform uniform(1, 10);
     dx_ = gaussian_pos.generate(rng);
     dy_ = gaussian_pos.generate(rng);
     // Log-normal prior for flux and bmaj
     logS_max_ = gaussian_logflux.generate(rng);
     logsize_ = gaussian_logsize.generate(rng);
	 lognu_max_ = uniform.generate(rng);
	 alpha_thick_ = gaussian_alpha_thick.generate(rng);
	 alpha_thin_ = gaussian_alpha_thin.generate(rng);
    // Dependent priors - from MOJAVE 15 GHz modelfits.
    //logflux_ = -1.74 + 1.62*rng.randn();
    //logbmaj_ = -0.71 - 0.44*logflux_ + 0.84*rng.randn();
    // Dependent priors - modified for 43 GHz.
//    logflux_ = -2.29 + 1.62*rng.randn();
//    logbmaj_ = -1.71 - 0.44*logflux_ + 0.84*rng.randn();
}

double JetGaussianComponent::perturb(DNest4::RNG &rng) {
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
		DNest4::Uniform uniform(1, 10);
		log_H += uniform.perturb(lognu_max_, rng);

	}
	else if(which == 5)
	{
		DNest4::Gaussian gaussian_alpha_thick(1.0, 0.5);
		log_H += gaussian_alpha_thick.perturb(alpha_thick_, rng);
	}
	else
	{
		DNest4::Gaussian gaussian_alpha_thin(-1.0, 0.5);
		log_H += gaussian_alpha_thin.perturb(alpha_thin_, rng);
	}
    return log_H;
}


CoreGaussianComponent::CoreGaussianComponent() {}

CoreGaussianComponent::CoreGaussianComponent(const CoreGaussianComponent &other)
{
	a_ = other.a_;
	PA_ = other.PA_;
	logsize_1_ = other.logsize_1_;
	k_r_ = other.k_r_;
	logS_1_ = other.logS_1_;
	alpha_ = other.alpha_;
}

CoreGaussianComponent *CoreGaussianComponent::clone()
{
	return new CoreGaussianComponent(*this);
}

double CoreGaussianComponent::get_logsize(double nu)
{
	return logsize_1_ - (1./k_r_)*log(nu);
}

double CoreGaussianComponent::get_logflux(double nu)
{
	return logS_1_ + alpha_*log(nu);
}

std::pair<double, double> CoreGaussianComponent::get_pos(double nu)
{
	double distance = a_*pow(nu, -1/k_r_);
	return {distance*sin(PA_), distance*cos(PA_)};
}

void CoreGaussianComponent::print(std::ostream &out) const
{
	out << a_ << "\t" << PA_ << "\t" << logsize_1_ << "\t" << k_r_ << "\t" << logS_1_ << "\t" << alpha_ << "\t";
}

std::string CoreGaussianComponent::description() const
{
	std::string descr;
	descr += "a PA logsize_1 k_r logS_1 alpha";
	return descr;
}

void CoreGaussianComponent::from_prior(DNest4::RNG &rng)
{
	DNest4::TruncatedCauchy cauchy_pos(0.0, 0.1, 0.0, 1.0);
	DNest4::Uniform gaussian_direction(0.0, M_PI);
	DNest4::Gaussian gaussian_logsize(-1.0, 2.0);
	DNest4::Gaussian gaussian_logflux(-1.0, 1.0);
	DNest4::Uniform uniform_alpha(-1.0, 1.0);
	DNest4::Uniform uniform_k(0.0, 3.0);
	a_ = cauchy_pos.generate(rng);
	PA_ = gaussian_direction.generate(rng);
	logsize_1_ = gaussian_logsize.generate(rng);
	k_r_ = uniform_k.generate(rng);
	logS_1_ = gaussian_logflux.generate(rng);
	alpha_ = uniform_alpha.generate(rng);
}

double CoreGaussianComponent::perturb(DNest4::RNG &rng)
{
	double log_H = 0.;
	int which = rng.rand_int(6);
	if(which == 0)
	{
		DNest4::TruncatedCauchy cauchy_pos(0.0, 0.1, 0.0, 1.0);
		log_H += cauchy_pos.perturb(a_, rng);
	}
	else if(which == 1)
	{
		DNest4::Uniform gaussian_direction(0.0, M_PI);
		log_H += gaussian_direction.perturb(PA_, rng);
	}
	else if(which == 2)
	{
		DNest4::Gaussian gaussian_logsize(-1.0, 2.0);
		log_H += gaussian_logsize.perturb(logsize_1_, rng);
	}
	else if(which == 3)
	{
		DNest4::Gaussian gaussian_logflux(-1.0, 1.0);
		log_H += gaussian_logflux.perturb(logS_1_, rng);
	}
	else if(which == 4)
	{
		DNest4::Uniform uniform_alpha(-1.0, 1.0);
		log_H += uniform_alpha.perturb(alpha_, rng);
	}
	else
	{
		DNest4::Uniform uniform_k(0.0, 3.0);
		log_H += uniform_k.perturb(k_r_, rng);
	}
	return 0;
}

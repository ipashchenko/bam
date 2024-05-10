#include "Utils.h"
#include "MyConditionalPrior.h"
#include "LogNormal.h"

using namespace DNest4;
extern const ComponentType component_type;


MyConditionalPrior::MyConditionalPrior()
{
	if(!raPrior){
		raPrior = std::make_shared<Gaussian>(0., 10.);
	}
	if(!decPrior){
		decPrior = std::make_shared<Gaussian>(0., 10.);
	}
	if(!FluxSizePrior){
		FluxSizePrior = std::make_shared<Gaussian2D>(-1.74, -0.71, -0.8, 1.62, 1.03);
	}
	if(!ePrior){
		ePrior = std::make_shared<DNest4::Kumaraswamy>(3.5, 2.5);
	}
	if(!bpaPrior){
		bpaPrior = std::make_shared<DNest4::Uniform>(0., M_PI);
	}
}

void MyConditionalPrior::from_prior(RNG& rng) {}

double MyConditionalPrior::perturb_hyperparameters(RNG& rng)
{
    double logH = 0.0;

//    int which = rng.rand_int(4);
//
//    if(which == 0)
//    {
//        const DNest4::Gaussian gauss1(-0.25, 0.05);
//        logH += gauss1.perturb(typical_flux, rng);
//    }
//    else if(which == 1)
//    {
//        const DNest4::Gaussian gauss2(1.0, 0.10);
//        logH += gauss2.perturb(dev_log_flux, rng);
//    }
//    else if(which == 2)
//    {
//        const DNest4::Gaussian gauss3(-1.0, 0.05);
//        logH += gauss3.perturb(typical_radius, rng);
//    }
//    else if(which == 3)
//    {
//        const DNest4::Gaussian gauss4(1.00, 0.10);
//        logH += gauss4.perturb(dev_log_radius, rng);
//    }
//	else if(which == 4)
//	{
//		const DNest4::Gaussian gauss5(3.5, 0.1);
//		logH += gauss5.perturb(typical_a, rng);
//	}
//	else if(which == 5)
//	{
//		const DNest4::Gaussian gauss6(2.5, 0.1);
//		logH += gauss6.perturb(typical_b, rng);
//	}

	
    return logH;
}


double MyConditionalPrior::log_pdf(const std::vector<double>& vec) const
{
    double logp = 0.0;

    logp += raPrior->log_pdf(vec[0]);
    logp += decPrior->log_pdf(vec[1]);
	
	std::vector<double> fluxsize(2);
	fluxsize[0] = log(vec[2]);
	fluxsize[1] = vec[3];
	logp += FluxSizePrior->log_pdf(fluxsize);
	
	if(component_type == elliptical)
	{
		logp += ePrior->log_pdf(vec[4]);
		logp += bpaPrior->log_pdf(vec[5]);
	}
	
    return logp;
}

void MyConditionalPrior::from_uniform(std::vector<double>& vec) const
{
    vec[0] = raPrior->cdf_inverse(vec[0]);
    vec[1] = decPrior->cdf_inverse(vec[1]);
	
	
	std::vector<double> u_fluxsize = {vec.begin() + 2, vec.end()};
	std::vector<double> fluxsize = FluxSizePrior->cdf_inverse(u_fluxsize);
	vec[2] = exp(fluxsize[0]);
	vec[3] = fluxsize[1];

	if(component_type == elliptical)
	{
		vec[4] = ePrior->cdf_inverse(vec[4]);
		vec[5] = bpaPrior->cdf_inverse(vec[5]);
	}

}


void MyConditionalPrior::to_uniform(std::vector<double>& vec) const
{
    vec[0] = raPrior->cdf(vec[0]);
    vec[1] = decPrior->cdf(vec[1]);
	
	std::vector<double> fluxsize(2);
	fluxsize[0] = log(vec[2]);
	fluxsize[1] = vec[3];
	std::vector<double> u_fluxsize = FluxSizePrior->cdf(fluxsize);
	vec[2] = u_fluxsize[0];
	vec[3] = u_fluxsize[1];
	
	if(component_type == elliptical)
	{
		vec[4] = ePrior->cdf(vec[4]);
		vec[5] = bpaPrior->cdf(vec[5]);
	}
	
}

void MyConditionalPrior::print(std::ostream& out) const
{
//    out << typical_flux << ' ' << dev_log_flux << ' ';
//    out << typical_radius << ' ' << dev_log_radius << ' ';
//	out << typical_a << ' ' << typical_b << ' ';
}

#include "MyConditionalPrior.h"
#include "LogNormal.h"


using namespace DNest4;

MyConditionalPrior::MyConditionalPrior()
{
	std::cout << "MyCondPrior ctor\n";
	// FIXME:
	if(!FluxSizePrior) {
		FluxSizePrior = std::make_shared<Gaussian2D>(-1.74, -0.71, -0.8, 1.62, 1.03);
	}
}

void MyConditionalPrior::from_prior(RNG& rng)
{
	std::cout << "MyCondPrior from_prior\n";
//	rho = -0.8;
//	sigma_1 = 1.62;
//	sigma_2 = 1.03;
//	mean = {-1.74, -0.71};
	
//    const DNest4::Gaussian gauss1(-0.25, 0.05);
//    typical_flux = gauss1.generate(rng);
//    const DNest4::Gaussian gauss2(1.00, 0.10);
//    dev_log_flux = gauss2.generate(rng);
//
//    const DNest4::Gaussian gauss3(-1.0, 0.05);
//    typical_radius = gauss3.generate(rng);
//    const DNest4::Gaussian gauss4(1.00, 0.10);
//    dev_log_radius = gauss4.generate(rng);
	
//	const DNest4::Gaussian gauss5(3.5, 0.1);
//	typical_a = gauss5.generate(rng);
//
//	const DNest4::Gaussian gauss6(2.5, 0.1);
//	typical_b = gauss6.generate(rng);
}

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

    // Position
	DNest4::Uniform uniform_x(x_min, x_max);
	DNest4::Uniform uniform_y(y_min, y_max);
    logp += uniform_x.log_pdf(vec[0]);
    logp += uniform_y.log_pdf(vec[1]);
	
	// x, y, flux, bmaj
	std::vector<double> fluxsize = {vec.begin() + 2, vec.end()};
	logp += FluxSizePrior->log_pdf(fluxsize);
//	Gaussian2D gaussian_2_d(mean, rho, sigma_1, sigma_2);
//	logp += gaussian_2_d.log_pdf({})
//
//
//    // Flux
//    // DNest4::Laplace laplace1(typical_flux, dev_log_flux);
//    LogNormal lognorm(typical_flux, dev_log_flux);
//    //logp += -log(vec[2]) + laplace1.log_pdf(log(vec[2]));
//    logp += lognorm.log_pdf(vec[2]);
//
//    // Radius
//    // DNest4::Laplace laplace2(typical_radius, dev_log_radius);
//    DNest4::Gaussian gauss2(typical_radius, dev_log_radius);
//    //logp += -log(vec[3]) + laplace2.log_pdf(log(vec[3]));
//    logp += gauss2.log_pdf(vec[3]);

//	// Elongation
//	DNest4::Kumaraswamy kumaraswamy(typical_a, typical_b);
//	logp += kumaraswamy.log_pdf(vec[4]);
	
    return logp;
}

void MyConditionalPrior::from_uniform(std::vector<double>& vec) const
{
//    std::cout << "to_uniform: begin " << vec[0] << " " << vec[1] << " " << vec[2] << " " << vec[3] << "\n";
    // Position
	DNest4::Uniform uniform_x(x_min, x_max);
	DNest4::Uniform uniform_y(y_min, y_max);
//    DNest4::Gaussian gaussx(0.0, std);
    vec[0] = uniform_x.cdf_inverse(vec[0]);
//    DNest4::Gaussian gaussy(0.0, std);
    vec[1] = uniform_y.cdf_inverse(vec[1]);
	
	
	std::vector<double> u_fluxsize = {vec.begin() + 2, vec.end()};
	std::vector<double> fluxsize = FluxSizePrior->cdf_inverse(u_fluxsize);
	vec[2] = fluxsize[0];
	vec[3] = fluxsize[1];
	
//    // Flux
//    // DNest4::Laplace laplace1(typical_flux, dev_log_flux);
//    LogNormal lognorm(typical_flux, dev_log_flux);
//    vec[2] = lognorm.cdf_inverse(vec[2]);
//
//    // Radius
//    // DNest4::Laplace laplace2(typical_radius, dev_log_radius);
//    DNest4::Gaussian gauss2(typical_radius, dev_log_radius);
//    vec[3] = gauss2.cdf_inverse(vec[3]);
	
//	// Elongation
//	DNest4::Kumaraswamy kumaraswamy(typical_a, typical_b);
//	vec[4] = kumaraswamy.cdf_inverse(vec[4]);
//
//	// BPA
//	DNest4::Uniform uniform(0., M_PI);
//	vec[5] = uniform.cdf_inverse(vec[5]);

}

void MyConditionalPrior::to_uniform(std::vector<double>& vec) const
{
    //std::cout << "to_uniform" << vec[0] << " " << vec[1] << " " << vec[2] << " " << vec[3] << "\n";
    // Position
	DNest4::Uniform uniform_x(x_min, x_max);
	DNest4::Uniform uniform_y(y_min, y_max);
//    DNest4::Gaussian gaussx(0.0, std);
    vec[0] = uniform_x.cdf(vec[0]);
//    DNest4::Gaussian gaussy(0.0, std);
    vec[1] = uniform_y.cdf(vec[1]);
	
	std::vector<double> fluxsize = {vec.begin() + 2, vec.end()};
	std::vector<double> u_fluxsize = FluxSizePrior->cdf(fluxsize);
	vec[2] = u_fluxsize[0];
	vec[3] = u_fluxsize[1];
	
	
//    // Flux
//    // DNest4::Laplace laplace1(typical_flux, dev_log_flux);
//    LogNormal lognorm(typical_flux, dev_log_flux);
//    vec[2] = lognorm.cdf(vec[2]);
//
//    // Radius
//    // DNest4::Laplace laplace2(typical_radius, dev_log_radius);
//    DNest4::Gaussian gauss2(typical_radius, dev_log_radius);
//    vec[3] = gauss2.cdf(vec[3]);
	
//	// Elongation
//	DNest4::Kumaraswamy kumaraswamy(typical_a, typical_b);
//	vec[4] = kumaraswamy.cdf(vec[4]);
//
//	// BPA
//	DNest4::Uniform uniform(0., M_PI);
//	vec[5] = uniform.cdf(vec[5]);
	
}

void MyConditionalPrior::print(std::ostream& out) const
{
//    out << typical_flux << ' ' << dev_log_flux << ' ';
//    out << typical_radius << ' ' << dev_log_radius << ' ';
//	out << typical_a << ' ' << typical_b << ' ';
}

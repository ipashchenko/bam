#include <iostream>
#include <Distributions/ContinuousDistribution.h>
#include "DNestModel.h"
#include "Data.h"


DNestModel::DNestModel()
{
	// Mapping from antenna numbers (ant_i/j) to position in vector of antennas.
	std::unordered_map<int, int>& antennas_map = Data::get_instance().get_antennas_map();
	const std::vector<int>& ant_i = Data::get_instance().get_ant_i();
	const std::vector<int>& ant_j = Data::get_instance().get_ant_j();
	for (int k=0; k<ant_i.size(); k++) {
		ant_ik.emplace_back(antennas_map[ant_i[k]]);
		ant_jk.emplace_back(antennas_map[ant_j[k]]);
	}
}


void DNestModel::setPriors()
{
	Jprior = make_prior<DNest4::Gaussian>(-5., 1.);
	Oprior = make_prior<DNest4::TruncatedCauchy>(1.0, 0.05, 0., 2.);
}

void DNestModel::from_prior(DNest4::RNG &rng) {
	
	setPriors();
	
	for(double & logj : per_antenna_logjitter)
	{
		logj = Jprior->generate(rng);
	}
	for(double & off : per_antenna_offset)
	{
		if(use_offsets)
		{
			off = Oprior->generate(rng);
		}
		else
		{
			off = 1.0;
		}
	}
	per_antenna_offset[reference_antenna_number] = 1.0;
	calculate_var();
	calculate_offset();
	// Assert #comp > 0
	do
	{	components.from_prior(rng);
	}while (components.get_components().size() == 0);
    components.consolidate_diff();
    calculate_sky_mu();
}


void DNestModel::calculate_var()
{
	auto antenna_map = Data::get_instance().get_antennas_map();
	const std::valarray<double>& sigma = Data::get_instance().get_sigma();
	// Measured variance
	std::valarray<double> measured_var = sigma*sigma;
	if(use_jitters)
	{
		// Jitter
		std::valarray<double> current_jitter(0.0, sigma.size());
		for (size_t k = 0; k < sigma.size(); k++)
		{
			current_jitter[k] = per_antenna_logjitter[antenna_map[ant_ik[k]]] + per_antenna_logjitter[antenna_map[ant_jk[k]]];
		}
		var = measured_var + exp(current_jitter);
	}
	else
	{
		var = measured_var;
	}
}


void DNestModel::calculate_offset()
{
	auto antenna_map = Data::get_instance().get_antennas_map();
	const std::valarray<double>& sigma = Data::get_instance().get_sigma();
	std::valarray<double> current_offset(0.0, sigma.size());
	for (size_t k=0; k<sigma.size(); k++) {
		current_offset[k] = per_antenna_offset[antenna_map[ant_ik[k]]] * per_antenna_offset[antenna_map[ant_jk[k]]];
	}
	offset = current_offset;
}

double DNestModel::perturb(DNest4::RNG &rng) {
    double logH = 0.;
	double u = rng.rand();
	
	double r_jitter = 0.0;
	double r_offset = 0.0;
	if (use_jitters) {
		r_jitter = 0.2;
	}
	if (use_offsets) {
		r_offset = 0.2;
	}
	
	
    // Perturb jitter
    if(u < r_jitter) {
	
		int which_antenna = rng.rand_int(per_antenna_logjitter.size());
		auto perturbed_jitter = per_antenna_logjitter[which_antenna];
		Jprior->perturb(perturbed_jitter, rng);
		per_antenna_logjitter[which_antenna] = perturbed_jitter;
        // No need to re-calculate model. Re-calculate variance and calculate loglike.
		calculate_var();
    }
	
	else if(r_jitter < u && u < (r_jitter + r_offset)) {
		int which_antenna = reference_antenna_number;
		while(which_antenna == reference_antenna_number)
		{
			which_antenna = rng.rand_int(per_antenna_offset.size());
		}
		auto perturbed_offset = per_antenna_offset[which_antenna];
		Oprior->perturb(perturbed_offset, rng);
		per_antenna_offset[which_antenna] = perturbed_offset;
		// No need to re-calculate model. Re-calculate variance and calculate loglike.
		calculate_offset();
	}
	
    // Perturb SkyModel
    else {
        logH += components.perturb(rng);
		// Do not allow zero components
		if(components.get_components().size() == 0)
		{
			return -std::numeric_limits<double>::max();
		}
		
        // Pre-reject
        if(rng.rand() >= exp(logH)) {
            return -1E300;
        }
        else
            logH = 0.0;
	
		components.consolidate_diff();
		// After this ``removed`` is empty and gone to ``added`` with ``-`` sign. We use ``added`` when
		// ``update = True``. Else we use all components.
        // This shouldn't be called in case of pre-rejection
        calculate_sky_mu();
    }
    return logH;
}


void DNestModel::calculate_sky_mu() {
    const std::valarray<double>& u = Data::get_instance().get_u();
    const std::valarray<double>& v = Data::get_instance().get_v();

    std::valarray<double> theta, pi_D_rho;
    double c;
    std::valarray<double> b;
    std::valarray<double> ft;

    // Update or from scratch?
	// From DNest4:
//	bool update_ = components.get_removed().size() == 0
    bool update = (components.get_added().size() < components.get_components().size()) &&
        (counter <= 20);
    // Get the components: comps will have only ``added`` if updating, all components if from scratch.
    const std::vector<std::vector<double>>& comps = (update)?(components.get_added()):
                                                    (components.get_components());

    if(!update)
    {
        // Zero the sky model prediction
        mu_real *= 0.0;
        mu_imag *= 0.0;
        //std::cout << "full " << counter << "\n";
        counter = 0;

    } else {
        //std::cout << "diff " << counter << "\n";
        counter++;
    }
	
	
	switch (component_type)
	{
	case circular:
		{
			for(const auto& comp: comps)
			{
				// Circular Gaussian
				// Phase shift due to not being in a phase center
				theta = 2*M_PI*mas_to_rad*(u*comp[0] + v*comp[1]);
				// Calculate FT of a Gaussian in a phase center
				c = pow(M_PI*exp(comp[3])*mas_to_rad, 2)/(4.*log(2.));
				ft = comp[2]*exp(-c*(u*u + v*v));
				// Prediction of visibilities
				mu_real += ft*cos(theta);
				mu_imag += ft*sin(theta);
			}
			break;
		}
		
	case sphere:
		{
			for(const auto& comp: comps)
			{
				// Circular Gaussian
				// Phase shift due to not being in a phase center
				theta = 2*M_PI*mas_to_rad*(u*comp[0] + v*comp[1]);
				// Calculate FT of a Sphere in a phase center
				pi_D_rho = M_PI*exp(comp[3])*mas_to_rad*sqrt(u*u + v*v);
				ft = 3.*comp[2]*(sin(pi_D_rho) - pi_D_rho*cos(pi_D_rho)) / pow(pi_D_rho, 3.);
				// Prediction of visibilities
				mu_real += ft*cos(theta);
				mu_imag += ft*sin(theta);
			}
			break;
		}
	
	case elliptical:
		{
			for(const auto& comp: comps)
			{
				// Elliptical Gaussian
				// Phase shift due to not being in a phase center
				theta = 2*M_PI*mas_to_rad*(u*comp[0] + v*comp[1]);
				// Calculate FT of a Gaussian in a phase center
				// NOTE: comp[3] is bmin!
				c = pow(M_PI*exp(comp[3])*mas_to_rad/comp[4], 2)/(4.*log(2.));
				b = comp[4]*comp[4]*pow((u*cos(comp[5]) - v*sin(comp[5])), 2) + pow((u*sin(comp[5]) + v*cos(comp[5])), 2);
				ft = comp[2]*exp(-c*b);
				// Prediction of visibilities
				mu_real += ft*cos(theta);
				mu_imag += ft*sin(theta);
			}
			break;
		}
	}
}


double DNestModel::log_likelihood() const {
    const std::valarray<double>& vis_real = Data::get_instance().get_vis_real();
    const std::valarray<double>& vis_imag = Data::get_instance().get_vis_imag();

    // Complex Gaussian sampling distribution
    std::valarray<double> result = -log(2*M_PI*var) - 0.5*(pow(vis_real - mu_real*offset, 2) + pow(vis_imag - mu_imag*offset, 2))/var;
    double loglik = result.sum();
    return loglik;
}


void DNestModel::print(std::ostream &out) const {
	if(use_jitters)
	{
		for (double logj : per_antenna_logjitter)
		{
			out << logj << '\t';
		}
	}
	if(use_offsets)
	{
		for (double off : per_antenna_offset)
		{
			out << off << '\t';
		}
	}
    components.print(out);
}


std::string DNestModel::description() const
{
	std::string descr;
	
	if (use_jitters) {
		for (int i = 0; i < per_antenna_logjitter.size(); i++) {
			descr += "logjitter[" + std::to_string(i) + "] ";
		}
	}
	
	if (use_offsets) {
		for (int i = 0; i < per_antenna_offset.size(); i++) {
			descr += "offset[" + std::to_string(i) + "] ";
		}
	}
	
	// The rest is all what happens when you call .print on an RJObject
	descr += " dim_components max_num_components ";
	
	// Then the hyperparameters (i.e. whatever MyConditionalPrior::print prints)
//    descr += " typical_flux dev_log_flux typical_radius dev_log_radius typical_a typical_b";
	
	// Then the actual number of components
	descr += " num_components ";
	
	// Then it's all the components, padded with zeros
	// max_num_components is known in this model, so that's how far the
	// zero padding goes.
	for (int i = 0; i < components.get_max_num_components(); ++i)
		descr += " x[" + std::to_string(i) + "] ";
	for (int i = 0; i < components.get_max_num_components(); ++i)
		descr += " y[" + std::to_string(i) + "] ";
	for (int i = 0; i < components.get_max_num_components(); ++i)
		descr += " flux[" + std::to_string(i) + "] ";
	for (int i = 0; i < components.get_max_num_components(); ++i)
		descr += " logsize[" + std::to_string(i) + "] ";
	
	if (component_type == elliptical)
	{
		for (int i = 0; i < components.get_max_num_components(); ++i)
			descr += " e[" + std::to_string(i) + "] ";
		for (int i = 0; i < components.get_max_num_components(); ++i)
			descr += " bpa[" + std::to_string(i) + "] ";
	}
    return descr;
}

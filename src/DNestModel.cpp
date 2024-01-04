#include <iostream>
#include <Distributions/ContinuousDistribution.h>
#include "DNestModel.h"
#include "Data.h"

// TODO: Optionally, use per-antenna jitter!
DNestModel::DNestModel() :
    counter(0),
	reference_antenna_number(4),
    components(6, 20, false, MyConditionalPrior(-10., 1., -7., 2.), DNest4::PriorType::log_uniform) {
	
	// Mapping from antenna numbers (ant_i/j) to position in vector of antennas.
	std::unordered_map<int, int>& antennas_map = Data::get_instance().get_antennas_map();
	const std::vector<int>& ant_i = Data::get_instance().get_ant_i();
	const std::vector<int>& ant_j = Data::get_instance().get_ant_j();
	for (int k=0; k<ant_i.size(); k++) {
		ant_ik.emplace_back(antennas_map[ant_i[k]]);
		ant_jk.emplace_back(antennas_map[ant_j[k]]);
	}
	
	int n_antennas = Data::get_instance().n_antennas();
	per_antenna_logjitter.resize(n_antennas);
	per_antenna_offset.resize(n_antennas);
}


void DNestModel::from_prior(DNest4::RNG &rng) {
    // I do this because in ``calculate_sky_mu`` ``mu_real`` and ``mu_imag`` are multiplied and added.
    const std::valarray<double>& u = Data::get_instance().get_u();
    std::valarray<double> zero (0.0, u.size());
    mu_real = zero;
    mu_imag = zero;
	for(double & logj : per_antenna_logjitter)
	{
		logj = -4.0 + 2.0*rng.randn();
	}
	DNest4::TruncatedCauchy truncated_cauchy(1.0, 0.05, 0., 2.);
	for(double & off : per_antenna_offset)
	{
		off = truncated_cauchy.generate(rng);
	}
	per_antenna_offset[reference_antenna_number] = 1.0;
	calculate_var();
	calculate_offset();
    components.from_prior(rng);
    components.consolidate_diff();
    calculate_sky_mu();
}


void DNestModel::calculate_var()
{
	auto antenna_map = Data::get_instance().get_antennas_map();
	const std::valarray<double>& sigma = Data::get_instance().get_sigma();
	// Measured variance
	std::valarray<double> measured_var = sigma*sigma;
	// Jitter
	std::valarray<double> current_jitter (0.0, sigma.size());
	for (size_t k=0; k<sigma.size(); k++) {
		current_jitter[k] = per_antenna_logjitter[antenna_map[ant_ik[k]]] + per_antenna_logjitter[antenna_map[ant_jk[k]]];
	}
	var = measured_var + exp(current_jitter);
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
	
    // Perturb jitter
    if(u <= 0.2) {
	
		int which_antenna = rng.rand_int(per_antenna_logjitter.size());
		auto perturbed_jitter = per_antenna_logjitter[which_antenna];
        logH -= -0.5*pow((perturbed_jitter + 4)/2.0, 2.0);
        perturbed_jitter += 2.0*rng.randh();
        logH += -0.5*pow((perturbed_jitter + 4)/2.0, 2.0);
		per_antenna_logjitter[which_antenna] = perturbed_jitter;
		

        // Pre-reject
        if(rng.rand() >= exp(logH)) {
            return -1E300;
        }
        else {
			logH = 0.0;
		}
        // No need to re-calculate model. Re-calculate variance and calculate loglike.
		calculate_var();
    }
	
	else if(0.2 <= u && u <= 0.4) {
		DNest4::TruncatedCauchy truncated_cauchy(1.0, 0.05, 0., 2.);
		int which_antenna = reference_antenna_number;
		while(which_antenna == reference_antenna_number)
		{
			which_antenna = rng.rand_int(per_antenna_offset.size());
		}
		auto perturbed_offset = per_antenna_offset[which_antenna];
		logH += truncated_cauchy.perturb(perturbed_offset, rng);
		per_antenna_offset[which_antenna] = perturbed_offset;
	
		// No need to re-calculate model. Re-calculate variance and calculate loglike.
		calculate_offset();
	}
	
    // Perturb SkyModel
    else {
        logH += components.perturb(rng);
        components.consolidate_diff();
        // After this ``removed`` is empty and gone to ``added`` with ``-`` sign. We use ``added`` when
        // ``update = True``. Else we use all components.

        // Pre-reject
        if(rng.rand() >= exp(logH)) {
            return -1E300;
        }
        else
            logH = 0.0;

        // This shouldn't be called in case of pre-rejection
        calculate_sky_mu();
    }
    return logH;
}


void DNestModel::calculate_sky_mu() {
    const std::valarray<double>& u = Data::get_instance().get_u();
    const std::valarray<double>& v = Data::get_instance().get_v();

    std::valarray<double> theta;
    double c;
    std::valarray<double> b;
    std::valarray<double> ft;

    // Update or from scratch?
    bool update = (components.get_added().size() < components.get_components().size()) &&
        (counter <= 20);
    // Get the components
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

    for(const auto& comp: comps)
    {
        // Phase shift due to not being in a phase center
        theta = 2*M_PI*mas_to_rad*(u*comp[0]+v*comp[1]);
        // Calculate FT of a Gaussian in a phase center
        c = pow(M_PI*exp(comp[3])*mas_to_rad, 2)/(4.*log(2.));
		b = comp[4]*comp[4]*pow((u*cos(comp[5]) - v*sin(comp[5])), 2) + pow((u*sin(comp[5]) + v*cos(comp[5])), 2);
        ft = comp[2] * exp(-c*b);
        // Prediction of visibilities
        mu_real += ft*cos(theta);
        mu_imag += ft*sin(theta);
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
	for(double logj : per_antenna_logjitter)
	{
		out << logj << '\t';
	}
	for(double off : per_antenna_offset)
	{
		out << off << '\t';
	}
    components.print(out);
}


std::string DNestModel::description() const
{
    std::string descr;

    // Anything printed by DNestModel::print (except the last line)
	for(int i = 0; i < per_antenna_logjitter.size(); i++)
	{
		descr += "logjitter[" + std::to_string(i) + "] ";
	}
	
	for(int i = 0; i < per_antenna_offset.size(); i++)
	{
		descr += "offset[" + std::to_string(i) + "] ";
	}
	
	
	// The rest is all what happens when you call .print on an RJObject
    descr += " dim_components max_num_components ";

    // Then the hyperparameters (i.e. whatever MyConditionalPrior::print prints)
    descr += " typical_flux dev_log_flux typical_radius dev_log_radius typical_a typical_b";

    // Then the actual number of components
    descr += " num_components ";

    // Then it's all the components, padded with zeros
    // max_num_components is known in this model, so that's how far the
    // zero padding goes.
    for(int i=0; i<components.get_max_num_components(); ++i)
        descr += " x[" + std::to_string(i) + "] ";
    for(int i=0; i<components.get_max_num_components(); ++i)
        descr += " y[" + std::to_string(i) + "] ";
    for(int i=0; i<components.get_max_num_components(); ++i)
        descr += " flux[" + std::to_string(i) + "] ";
    for(int i=0; i<components.get_max_num_components(); ++i)
        descr += " logbmaj[" + std::to_string(i) + "] ";
	for(int i=0; i<components.get_max_num_components(); ++i)
		descr += " e[" + std::to_string(i) + "] ";
	for(int i=0; i<components.get_max_num_components(); ++i)
		descr += " bpa[" + std::to_string(i) + "] ";
    return descr;
}

#include <iostream>
#include <Distributions/ContinuousDistribution.h>
#include "DNestModel.h"
#include "Data.h"
#include "Component.h"


DNestModel::DNestModel() : use_logjitter(true), use_speedup(true), component_ft_counter(0)
{
	size_t number_of_jet_components = 6;
    sky_model = new SkyModel(number_of_jet_components);
	old_sky_model = sky_model->clone();
}


DNestModel::~DNestModel()
{
    delete sky_model;
	delete old_sky_model;
}


DNestModel::DNestModel(const DNestModel& other)
{
    sky_model = new SkyModel(*other.sky_model);
	old_sky_model = new SkyModel(*other.old_sky_model);
    logjitter = other.logjitter;
    use_logjitter = other.use_logjitter;
	use_speedup = other.use_speedup;
	sky_model_mu = other.sky_model_mu;
	mu_full = other.mu_full;
	jet_origin_x = other.jet_origin_x;
	jet_origin_y = other.jet_origin_y;
	component_ft_counter = other.component_ft_counter;
}


DNestModel& DNestModel::operator=(const DNestModel& other)
	{
    if (this != &other)
	{
        *(sky_model) = *(other.sky_model);
        *(old_sky_model) = *(other.old_sky_model);
        logjitter = other.logjitter;
        use_logjitter = other.use_logjitter;
		use_speedup = other.use_speedup;
		sky_model_mu = other.sky_model_mu;
		mu_full = other.mu_full;
		jet_origin_x = other.jet_origin_x;
		jet_origin_y = other.jet_origin_y;
		component_ft_counter = other.component_ft_counter;
    }
    return *this;
}


void DNestModel::from_prior(DNest4::RNG &rng)
{
	DNest4::Cauchy cauchy_origin(0.0, 0.1);
	DNest4::Gaussian gaussian_jitter(-4.0, 2.0);
	const std::unordered_map<std::string, double> band_freq_map = Data::get_instance().get_band_freq_map();
	for (const auto& [band, freq] : band_freq_map)
	{
		// I do this because in ``calculate_sky_mu`` ``sky_model_mu_real`` and ``sky_model_mu_imag`` are multiplied and added.
		const ArrayXd &u = Data::get_instance().get_u(band);
		ArrayXcd zero = ArrayXcd::Zero(u.size());
		sky_model_mu[band] = zero;
		mu_full[band] = zero;
		jet_origin_x[band] = cauchy_origin.generate(rng);
		jet_origin_y[band] = cauchy_origin.generate(rng);
		logjitter[band] = gaussian_jitter.generate(rng);
		
	}
    sky_model->from_prior(rng);
    calculate_sky_mu(false);
	// old SkyModel now is the copy of the current one
	old_sky_model = sky_model->clone();
	shift_sky_mu();
}


double DNestModel::perturb(DNest4::RNG &rng)
{
    double logH = 0.;

    double u = rng.rand();
    double r_logjitter = 0.0;
    if (use_logjitter)
	{
        r_logjitter = 0.1;
    }

    // Perturb jitter
    if(u <= r_logjitter)
	{
		
		std::vector<std::string> bands = Data::get_instance().get_bands();
		int n_bands = bands.size();
		int which = rng.rand_int(n_bands);
		std::string band = bands[which];
		double logjitter_in_band;
		DNest4::Gaussian gaussian_jitter(-4.0, 2.0);
		logjitter_in_band = logjitter[band];
		logH += gaussian_jitter.perturb(logjitter_in_band, rng);
		logjitter[band] = logjitter_in_band;
        // No need to re-calculate model. Just calculate loglike.
    }

    // Perturb SkyModel
    else if(u >= r_logjitter && u < r_logjitter + 0.7)
	{
        logH += sky_model->perturb(rng);
		// Now the old sky_model has the same boolean perturbed vector
		old_sky_model->set_perturbed(sky_model->get_perturbed());
        // Pre-reject
        if(rng.rand() >= exp(logH))
		{
            return -1E300;
        }
        else
            logH = 0.0;

        // This shouldn't be called in case of pre-rejection
        calculate_sky_mu(true);
		shift_sky_mu();
    }
	
	// Perturb per-band phase centers
	else
	{
		std::vector<std::string> bands = Data::get_instance().get_bands();
		int n_bands = bands.size();
		int which = rng.rand_int(n_bands);
		std::string band = bands[which];
		DNest4::Cauchy cauchy_origin(0.0, 0.1);
		int which_xy = rng.rand_int(2);
		double origin;
		if(which_xy == 0)
		{
			origin = jet_origin_x[band];
			logH += cauchy_origin.perturb(origin, rng);
			jet_origin_x[band] = origin;
		}
		else
		{
			origin = jet_origin_y[band];
			logH += cauchy_origin.perturb(origin, rng);
			jet_origin_y[band] = origin;
		}
		shift_sky_mu();
	}
    return logH;
}


void DNestModel::calculate_sky_mu(bool update)
{
	const std::unordered_map<std::string, double> band_freq_map = Data::get_instance().get_band_freq_map();
	for (const auto& [band, freq] : band_freq_map)
	{
		const ArrayXd &u = Data::get_instance().get_u(band);
		const ArrayXd &v = Data::get_instance().get_v(band);
		
		if(use_speedup && update && component_ft_counter < 30)
		{
//			std::cout << "Speed up with counter = " << component_ft_counter << "\n";
			// Switching boolean perturbed[i] to false here
			sky_model->ft_from_perturbed(freq, u, v);
			//! Here all components in old_sky_model are not perturbed!
			old_sky_model->ft_from_perturbed(freq, u, v);
			sky_model_mu[band] = sky_model->get_mu() - old_sky_model->get_mu();
			component_ft_counter += 1;
		}
		else
		{
//			std::cout << "Full : resetting counter!\n";
			sky_model->ft_from_all(freq, u, v);
			component_ft_counter = 0;
			sky_model_mu[band] = sky_model->get_mu();
		}
	}
}

//* I need to decouple original predictions of SkyModel and shifted predictions. Because I change it each perturb!
void DNestModel::shift_sky_mu()
{
	const std::unordered_map<std::string, double> band_freq_map = Data::get_instance().get_band_freq_map();
	const std::complex<double> j(0.0, 1.0);
	for (const auto& [band, freq] : band_freq_map)
	{
		ArrayXd &u = Data::get_instance().get_u(band);
		ArrayXd &v = Data::get_instance().get_v(band);
		ArrayXcd theta = 2 * M_PI * j * mas_to_rad * (u * jet_origin_x[band] + v * jet_origin_y[band]);
		mu_full[band] = sky_model_mu[band] * exp(theta);
	}
}

double DNestModel::log_likelihood()
{
	double loglik = 0.0;
	const std::unordered_map<std::string, double> band_freq_map = Data::get_instance().get_band_freq_map();
	for (const auto& [band, freq] : band_freq_map)
	{
			const ArrayXcd &vis_obs = Data::get_instance().get_vis_real(band);
			const ArrayXd &sigma = Data::get_instance().get_sigma(band);
			
			ArrayXd var = sigma*sigma;
			
			if (use_logjitter)
			{
				var = var + exp(2.0*logjitter[band]);
			}
			
			// Complex Gaussian sampling distribution
			ArrayXd result = -log(2*M_PI*var) - 0.5*square(abs(vis_obs - mu_full[band]))/var;
			loglik += result.sum();
	}
    return loglik;
}


void DNestModel::print(std::ostream &out) const
{
    if(use_logjitter)
	{
		for (const auto& [band, logjitter_band] : logjitter)
		{
			out << logjitter_band << "\t";
		}
    }
	for (const auto& [band, origin] : jet_origin_x)
	{
		out << origin << "\t";
	}
	for (const auto& [band, origin] : jet_origin_y)
	{
		out << origin << "\t";
	}
	sky_model->print(out);
}


std::string DNestModel::description() const
{
    std::string descr;

    if(use_logjitter)
	{
		for (const auto& [band, logjitter_band] : logjitter)
		{
			descr += "logjitter_" + band + "\t";
		}
    }
	for (const auto& [band, origin] : jet_origin_x)
	{
		descr += "x_origin_" + band + "\t";
	}
	for (const auto& [band, origin] : jet_origin_y)
	{
		descr += "y_origin_" + band + "\t";
	}
	descr += sky_model->description();
    return descr;
}

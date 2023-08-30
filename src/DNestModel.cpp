#include <iostream>
#include <Distributions/ContinuousDistribution.h>
#include "DNestModel.h"
#include "Data.h"
#include "Component.h"

#ifdef NDEBUG
#define DEBUG(x)
#else
#define DEBUG(x) do { std::cerr << x << std::endl; } while (0)
#endif


DNestModel::DNestModel() : use_logjitter(true), use_speedup(true), use_image_shifts_from_core(true), component_ft_counter(0)
{
	size_t number_of_jet_components = 2;
    sky_model = new SkyModel(number_of_jet_components);
	old_sky_model = sky_model->clone();
	DEBUG("DM ctor sky_model : " + sky_model->print());
	DEBUG("DM ctor old_sky_model : " + old_sky_model->print());
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
	use_image_shifts_from_core = other.use_image_shifts_from_core;
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
		use_image_shifts_from_core = other.use_image_shifts_from_core;
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
	DNest4::Gaussian gaussian_origin(0.0, 0.1);
	DNest4::Gaussian gaussian_jitter(-4.0, 2.0);
	const std::unordered_map<std::string, double> band_freq_map = Data::get_instance().get_band_freq_map();
	for (const auto& [band, freq] : band_freq_map)
	{
		// I do this because in ``calculate_sky_mu`` ``sky_model_mu_real`` and ``sky_model_mu_imag`` are multiplied and added.
		const ArrayXd &u = Data::get_instance().get_u(band);
		ArrayXcd zero = ArrayXcd::Zero(u.size());
		sky_model_mu[band] = zero;
		mu_full[band] = zero;
		jet_origin_x[band] = gaussian_origin.generate(rng);
		jet_origin_y[band] = gaussian_origin.generate(rng);
		logjitter[band] = gaussian_jitter.generate(rng);
		
	}
    sky_model->from_prior(rng);
	// old SkyModel now is the copy of the current one
	old_sky_model = sky_model->clone();
	DEBUG("DM from_prior sky_model : " + sky_model->print());
	DEBUG("DM from_prior old_sky_model : " + old_sky_model->print());
	calculate_sky_mu(false);
	shift_sky_mu();
}


double DNestModel::perturb(DNest4::RNG &rng)
{
    double logH = 0.;

    double u = rng.rand();
    double r_logjitter = 0.0;
    if (use_logjitter)
	{
        r_logjitter = 0.2;
    }

    // Perturb jitter
    if(u <= r_logjitter)
	{
		DEBUG("Perturbing Jitter");
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
    else if(u >= r_logjitter && u < r_logjitter + 0.5)
	{
		DEBUG("Perturbing SkyModel");
		DEBUG("Perturb: sky_model before perturb : " + sky_model->print());
		logH += sky_model->perturb(rng);
		DEBUG("Perturb: sky_model after perturb : " + sky_model->print());
		DEBUG("Perturb: old_sky_model after perturb : " + old_sky_model->print());
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
		
		// Reset perturbed flag here
		sky_model->reset_perturbed();
		old_sky_model->reset_perturbed();
		DEBUG("DNestModel.calculate_sky_mu - reset both sky_model counters");
		
		// Some checks
		auto old_perturbed = old_sky_model->get_perturbed();
		auto perturbed = sky_model->get_perturbed();
		if (std::find(std::begin(perturbed), std::end(perturbed), true) == std::end(perturbed)) // All false
		{
			DEBUG("All perturbed = false after resetting in sky_model");
		}
		else
		{
			DEBUG("All perturbed != false after resetting in sky_model!!!");
		}
		if (std::find(std::begin(old_perturbed), std::end(old_perturbed), true) == std::end(old_perturbed)) // All false
		{
			DEBUG("All perturbed = false after resetting in old_sky_model");
		}
		else
		{
			DEBUG("All perturbed != false after resetting in old_sky_model");
		}
		
		// Make old - copy of the current
		old_sky_model = sky_model->clone();
		DEBUG("Perturb: sky_model after cloning : " + sky_model->print());
		DEBUG("Perturb: old_sky_model after cloning : " + old_sky_model->print());
    }
	
	// Perturb per-band phase centers
	else
	{
		std::vector<std::string> bands = Data::get_instance().get_bands();
		int n_bands = bands.size();
		int which = rng.rand_int(n_bands);
		std::string band = bands[which];
		DEBUG("Perturbing phase centers band = " + band);
		DNest4::Gaussian gaussian_origin(0.0, 0.1);
		double origin;

		double uu = rng.rand();
		// Perturb both coordinates
		if(uu > 0.5)
		{
			origin = jet_origin_x[band];
			logH += gaussian_origin.perturb(origin, rng);
			jet_origin_x[band] = origin;
			origin = jet_origin_y[band];
			logH += gaussian_origin.perturb(origin, rng);
			jet_origin_y[band] = origin;
		}
		// Perturb only one coordinate
		else
		{
			int which_xy = rng.rand_int(2);
			if(which_xy == 0)
			{
				origin = jet_origin_x[band];
				logH += gaussian_origin.perturb(origin, rng);
				jet_origin_x[band] = origin;
			}
			else
			{
				origin = jet_origin_y[band];
				logH += gaussian_origin.perturb(origin, rng);
				jet_origin_y[band] = origin;
			}
		}
		shift_sky_mu();
	}
    return logH;
}


void DNestModel::calculate_sky_mu(bool update)
{
	DEBUG("Starting DM:calculate_sky_mu with ft_counter = " + std::to_string(component_ft_counter));
	const std::unordered_map<std::string, double> band_freq_map = Data::get_instance().get_band_freq_map();
	for (const auto& [band, freq] : band_freq_map)
	{
		const ArrayXd &u = Data::get_instance().get_u(band);
		const ArrayXd &v = Data::get_instance().get_v(band);
		
		if(use_speedup && update && component_ft_counter < 30)
		{
			DEBUG("DNestModel.calculate_sky_mu - Speed up with counter = " + std::to_string(component_ft_counter));
			// Change the model prediction by the difference between predictions of the new ond old perturbed component
			// FIXME: ???
			sky_model_mu[band] += (sky_model->ft_from_perturbed(freq, u, v) - old_sky_model->ft_from_perturbed(freq, u, v)).eval();
		}
		else
		{
			DEBUG("DNestModel.calculate_sky_mu - Full : will reset counter!");
			sky_model_mu[band] = sky_model->ft_from_all(freq, u, v);
		}
	}
	
	if(use_speedup && update && component_ft_counter < 30)
	{
		component_ft_counter += 1;
	}
	else
	{
		component_ft_counter = 0;
	}
}

//* I need to decouple original predictions of SkyModel and shifted predictions. Because I change it each perturb!
// TODO: Shift only if phase center of given band has changed or if sky_model has changed.
void DNestModel::shift_sky_mu()
{
	DEBUG("Shifting sky mu");
	const std::unordered_map<std::string, double> band_freq_map = Data::get_instance().get_band_freq_map();
	const std::complex<double> j(0.0, 1.0);
	for (const auto& [band, freq] : band_freq_map)
	{
		std::pair<double, double> reference{0., 0.};
		// Obtain core position at given band relative to jet origin
		if(use_image_shifts_from_core){
			reference = sky_model->get_core_position(freq);
		}
		ArrayXd &u = Data::get_instance().get_u(band);
		ArrayXd &v = Data::get_instance().get_v(band);
		mu_full[band] = sky_model_mu[band] * exp(2 * M_PI * j * mas_to_rad * (u*(-reference.first + jet_origin_x[band]) +
																			  v*(-reference.second + jet_origin_y[band])));
	}
}

double DNestModel::log_likelihood()
{
	DEBUG("Calculating logL");
	double loglik = 0.0;
	const std::unordered_map<std::string, double> band_freq_map = Data::get_instance().get_band_freq_map();
	for (const auto& [band, freq] : band_freq_map)
	{
		const ArrayXcd &vis_obs = Data::get_instance().get_vis_real(band);
		const ArrayXd &sigma = Data::get_instance().get_sigma(band);
		
		ArrayXd var = sigma*sigma;
		
		if (use_logjitter)
		{
			var = (var + exp(2.0*logjitter[band])).eval();
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

#include <iostream>
#include <Distributions/ContinuousDistribution.h>
#include "DNestModel.h"
#include "Data.h"
#include "Component.h"


DNestModel::DNestModel() : use_logjitter(true)
{
    sky_model = new SkyModel();
    int n_jetcomp = 5;
	auto* comp = new CoreGaussianComponent();
	sky_model->add_component(comp);
    for (int i=0; i<n_jetcomp; i++) {
        auto* comp = new JetGaussianComponent();
        sky_model->add_component(comp);
    }
}


DNestModel::~DNestModel()
{
    delete sky_model;
}


DNestModel::DNestModel(const DNestModel& other)
{
    sky_model = new SkyModel(*other.sky_model);
    logjitter = other.logjitter;
    use_logjitter = other.use_logjitter;
	sky_model_mu_real = other.sky_model_mu_real;
	sky_model_mu_imag = other.sky_model_mu_imag;
	mu_real_full = other.mu_real_full;
	mu_imag_full = other.mu_imag_full;
	jet_origin_x = other.jet_origin_x;
	jet_origin_y = other.jet_origin_y;
}


DNestModel& DNestModel::operator=(const DNestModel& other)
	{
    if (this != &other)
	{
        *(sky_model) = *(other.sky_model);
        logjitter = other.logjitter;
        use_logjitter = other.use_logjitter;
		sky_model_mu_real = other.sky_model_mu_real;
		sky_model_mu_imag = other.sky_model_mu_imag;
		mu_real_full = other.mu_real_full;
		mu_imag_full = other.mu_imag_full;
		jet_origin_x = other.jet_origin_x;
		jet_origin_y = other.jet_origin_y;
    }
    return *this;
}


// TODO: Here initialize containers with old and current per-component predictions.
void DNestModel::from_prior(DNest4::RNG &rng)
{
	DNest4::Gaussian gaussian_origin(0.0, 0.25);
	DNest4::Gaussian gaussian_jitter(-4.0, 2.0);
	const std::unordered_map<std::string, double> band_freq_map = Data::get_instance().get_band_freq_map();
	for (const auto& [band, freq] : band_freq_map)
	{
		// I do this because in ``calculate_sky_mu`` ``sky_model_mu_real`` and ``sky_model_mu_imag`` are multiplied and added.
		const std::valarray<double> &u = Data::get_instance().get_u(band);
		std::valarray<double> zero(0.0, u.size());
		sky_model_mu_real[band] = zero;
		sky_model_mu_imag[band] = zero;
		mu_real_full[band] = zero;
		mu_imag_full[band] = zero;
		jet_origin_x[band] = gaussian_origin.generate(rng);
		jet_origin_y[band] = gaussian_origin.generate(rng);
		logjitter[band] = gaussian_jitter.generate(rng);
		
	}
    sky_model->from_prior(rng);
    calculate_sky_mu();
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

        // Pre-reject
        if(rng.rand() >= exp(logH))
		{
            return -1E300;
        }
        else
            logH = 0.0;

        // This shouldn't be called in case of pre-rejection
        calculate_sky_mu();
		shift_sky_mu();
    }
	
	// Perturb per-band phase centers
	else
	{
		std::vector<std::string> bands = Data::get_instance().get_bands();
		int n_bands = bands.size();
		int which = rng.rand_int(n_bands);
		std::string band = bands[which];
		DNest4::Gaussian gaussian_origin(0.0, 0.25);
		int which_xy = rng.rand_int(2);
		double origin;
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
		shift_sky_mu();
	}
    return logH;
}


void DNestModel::calculate_sky_mu()
{
	const std::unordered_map<std::string, double> band_freq_map = Data::get_instance().get_band_freq_map();
	for (const auto& [band, freq] : band_freq_map)
	{
		const std::valarray<double> &u = Data::get_instance().get_u(band);
		const std::valarray<double> &v = Data::get_instance().get_v(band);
		
		// FT (calculate SkyModel prediction)
		sky_model->ft_from_all(freq, u, v);
		sky_model_mu_real[band] = sky_model->get_mu_real();
		sky_model_mu_imag[band] = sky_model->get_mu_imag();
	}
}

// FIXME: I need de-couple original predictions of SkyModel and shifted predictions. Because I change it each perturb!
void DNestModel::shift_sky_mu()
{
	const std::unordered_map<std::string, double> band_freq_map = Data::get_instance().get_band_freq_map();
	for (const auto& [band, freq] : band_freq_map)
	{
		const std::valarray<double> &u = Data::get_instance().get_u(band);
		const std::valarray<double> &v = Data::get_instance().get_v(band);
		std::valarray<double> theta = 2 * M_PI * mas_to_rad * (u * jet_origin_x[band] + v * jet_origin_y[band]);
		mu_real_full[band] = cos(theta) * sky_model_mu_real[band] - sin(theta) * sky_model_mu_imag[band];
		mu_imag_full[band] = cos(theta) * sky_model_mu_imag[band] + sin(theta) * sky_model_mu_real[band];
	}
}

double DNestModel::log_likelihood()
{
	double loglik = 0.0;
	const std::unordered_map<std::string, double> band_freq_map = Data::get_instance().get_band_freq_map();
	for (const auto& [band, freq] : band_freq_map)
	{
			const std::valarray<double> &vis_real = Data::get_instance().get_vis_real(band);
			const std::valarray<double> &vis_imag = Data::get_instance().get_vis_imag(band);
			const std::valarray<double> &sigma = Data::get_instance().get_sigma(band);
			
			std::valarray<double> var = sigma*sigma;
			
			if (use_logjitter)
			{
				var = var + exp(2.0*logjitter[band]);
			}
			
			// Complex Gaussian sampling distribution
			std::valarray<double> result = -log(2*M_PI*var) - 0.5*(pow(vis_real - mu_real_full[band], 2) +
				pow(vis_imag - mu_imag_full[band], 2))/var;
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
    sky_model->print(out);
	for (const auto& [band, origin] : jet_origin_x)
	{
		out << origin << "\t";
	}
	for (const auto& [band, origin] : jet_origin_y)
	{
		out << origin << "\t";
	}
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

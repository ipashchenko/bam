#include <iostream>
#include <Distributions/ContinuousDistribution.h>
#include "DNestModel.h"
#include "Data.h"
#include "Component.h"


DNestModel::DNestModel() : logjitter(0.0), use_logjitter(true)
{
    sky_model = new SkyModel();
    int n_jetcomp = 2;
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
    mu_real = other.mu_real;
    mu_imag = other.mu_imag;
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
        mu_real = other.mu_real;
        mu_imag = other.mu_imag;
		jet_origin_x = other.jet_origin_x;
		jet_origin_y = other.jet_origin_y;
    }
    return *this;
}


// TODO: Here initialize containers with old and current per-component predictions.
void DNestModel::from_prior(DNest4::RNG &rng)
{
	DNest4::Gaussian gaussian_origin(0.0, 0.25);
	const std::unordered_map<std::string, double> band_freq_map = Data::get_instance().get_band_freq_map();
	for (const auto& [band, freq] : band_freq_map)
	{
		// I do this because in ``calculate_sky_mu`` ``mu_real`` and ``mu_imag`` are multiplied and added.
		const std::valarray<double> &u = Data::get_instance().get_u(band);
		std::valarray<double> zero(0.0, u.size());
		mu_real[band] = zero;
		mu_imag[band] = zero;
		jet_origin_x[band] = gaussian_origin.generate(rng);
		jet_origin_y[band] = gaussian_origin.generate(rng);
	}
    logjitter = -4.0 + 2.0*rng.randn();
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
        logH -= -0.5*pow((logjitter+4)/2.0, 2.0);
        logjitter += 2.0*rng.randh();
        logH += -0.5*pow((logjitter+4)/2.0, 2.0);

        // Pre-reject
        if(rng.rand() >= exp(logH))
		{
            return -1E300;
        }
        else
            logH = 0.0;
        // No need to re-calculate model. Just calculate loglike.
    }

    // Perturb SkyModel
    else if(u >= r_logjitter && u < r_logjitter + 0.8)
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
		mu_real[band] = sky_model->get_mu_real();
		mu_imag[band] = sky_model->get_mu_imag();
	}
}

void DNestModel::shift_sky_mu()
{
	const std::unordered_map<std::string, double> band_freq_map = Data::get_instance().get_band_freq_map();
	for (const auto& [band, freq] : band_freq_map)
	{
		const std::valarray<double> &u = Data::get_instance().get_u(band);
		const std::valarray<double> &v = Data::get_instance().get_v(band);
		std::valarray<double> theta = 2 * M_PI * mas_to_rad * (u * jet_origin_x[band] + v * jet_origin_y[band]);
		mu_real[band] = cos(theta) * mu_real[band] - sin(theta) * mu_imag[band];
		mu_imag[band] = cos(theta) * mu_imag[band] + sin(theta) * mu_real[band];
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
			
			std::valarray<double> var = sigma * sigma;
			
			if (use_logjitter)
			{
				var = var + exp(2.0 * logjitter);
			}
			
			// Complex Gaussian sampling distribution
			std::valarray<double> result = -log(2 * M_PI * var) - 0.5 * (pow(vis_real - mu_real[band], 2) + pow(vis_imag - mu_imag[band], 2)) / var;
			loglik += result.sum();
	}
    return loglik;
}


void DNestModel::print(std::ostream &out) const
{
    if(use_logjitter)
	{
        out << logjitter << '\t';
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

    // Anything printed by DNestModel::print (except the last line)
    if(use_logjitter)
	{
        descr += "logjitter\t";
    }
    descr += sky_model->description();
	for (const auto& [band, origin] : jet_origin_x)
	{
		descr += "x_origin_" + band + "\t";
	}
	for (const auto& [band, origin] : jet_origin_y)
	{
		descr += "y_origin_" + band + "\t";
	}
	descr.pop_back();
    return descr;
}

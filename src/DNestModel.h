#ifndef BAM_DNESTMODEL_H
#define BAM_DNESTMODEL_H

#include <map>
#include "SkyModel.h"
#include "RNG.h"
#include "DNest4.h"


class DNestModel {
    public:

        DNestModel();
        ~DNestModel();
        DNestModel(const DNestModel& other);
        DNestModel& operator=(const DNestModel& other);

        // Generate the point from the prior
        void from_prior(DNest4::RNG& rng);

        // Metropolis-Hastings proposals
        double perturb(DNest4::RNG& rng);

        // Likelihood function
        double log_likelihood();

        // Print to stream
        void print(std::ostream& out) const;

        // Return string with column information
        std::string description() const;

    private:
		std::map<std::string, double> logjitter;
        bool use_logjitter;
		// Per-band position of the jet base relative to the phase center
		std::map<std::string, double> jet_origin_x;
		std::map<std::string, double> jet_origin_y;
        SkyModel* sky_model;
        SkyModel* old_sky_model;
        std::vector<bool> perturbed_components;
        // Prediction of SkyModel only
        std::unordered_map<std::string, ArrayXcd> sky_model_mu;
		// Predictions with shift
		std::unordered_map<std::string, ArrayXcd> mu_full;
        // This runs ``ft`` method of SkyModel with (u, v) from Data and updates SkyModel predictions
        void calculate_sky_mu(bool update);
		void shift_sky_mu();
		size_t component_ft_counter;
		bool use_speedup;
};

#endif //BAM_DNESTMODEL_H

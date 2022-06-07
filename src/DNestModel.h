#ifndef BAM_DNESTMODEL_H
#define BAM_DNESTMODEL_H

#include <valarray>
#include "SkyModel.h"
#include "RNG.h"
#include "DNest4.h"
#include "MyConditionalPrior.h"


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
        double log_likelihood() const;

        // Print to stream
        void print(std::ostream& out) const;

        // Return string with column information
        std::string description() const;

    private:
        double logjitter;
        bool use_logjitter;
        bool use_speedup;
        size_t ft_calc_counter;
        SkyModel* sky_model{};
        // Prediction of SkyModel only
        std::valarray<double> mu_real;
        std::valarray<double> mu_imag;
        // This runs ``ft`` method of SkyModel with (u, v) from Data and updates SkyModel predictions
        void calculate_sky_mu(bool update);
};

#endif //BAM_DNESTMODEL_H

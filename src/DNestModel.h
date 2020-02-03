#ifndef BAM_DNESTMODEL_H
#define BAM_DNESTMODEL_H

#include <valarray>
#include "RNG.h"
#include "DNest4.h"
#include "MyConditionalPrior.h"

const double mas_to_rad = 4.84813681109536e-09;

class DNestModel {
    public:

        DNestModel();

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
        unsigned int counter;
        DNest4::RJObject<MyConditionalPrior> components;
        std::valarray<double> mu_real;
        std::valarray<double> mu_imag;
        void calculate_sky_mu();
};

#endif //BAM_DNESTMODEL_H

#ifndef BAM__MYCONDITIONALPRIOR_H_
#define BAM__MYCONDITIONALPRIOR_H_

#pragma once

#include "DNest4.h"
#include "Gaussian2D.h"


// Hyperparameters setting interim prior for components properties
class MyConditionalPrior:public DNest4::ConditionalPrior
{
    private:

        // Parameters of hyper-distributions
        double x_min{-1.}, x_max{15.}, y_min{-3.}, y_max{3.};
		
        double perturb_hyperparameters(DNest4::RNG& rng);

    public:
        MyConditionalPrior();
		
		std::shared_ptr<Gaussian2D> FluxSizePrior;
		
        void from_prior(DNest4::RNG& rng);

        double log_pdf(const std::vector<double>& vec) const;
        void from_uniform(std::vector<double>& vec) const;
        void to_uniform(std::vector<double>& vec) const;

        void print(std::ostream& out) const;
        static const int weight_parameter = 2;

};

#endif //BAM__MYCONDITIONALPRIOR_H_

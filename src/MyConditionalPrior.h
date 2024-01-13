#ifndef BAM__MYCONDITIONALPRIOR_H_
#define BAM__MYCONDITIONALPRIOR_H_

#pragma once

#include "DNest4.h"
#include "Gaussian2D.h"


// Hyperparameters setting interim prior for components properties
class MyConditionalPrior:public DNest4::ConditionalPrior
{
    private:

        double perturb_hyperparameters(DNest4::RNG& rng);

    public:
        MyConditionalPrior();
		
		std::shared_ptr<DNest4::Uniform> raPrior;
		std::shared_ptr<DNest4::Uniform> decPrior;
		std::shared_ptr<Gaussian2D> FluxSizePrior;
		std::shared_ptr<DNest4::Kumaraswamy> ePrior;
		std::shared_ptr<DNest4::Uniform> bpaPrior;
		
        void from_prior(DNest4::RNG& rng);

        double log_pdf(const std::vector<double>& vec) const;
        void from_uniform(std::vector<double>& vec) const;
        void to_uniform(std::vector<double>& vec) const;

        void print(std::ostream& out) const;
        static const int weight_parameter = 2;

};

#endif //BAM__MYCONDITIONALPRIOR_H_

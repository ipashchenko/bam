#ifndef BAM_DNESTMODEL_H
#define BAM_DNESTMODEL_H

#pragma once


#include <valarray>
#include "RNG.h"
#include "DNest4.h"
#include "Utils.h"
#include "Data.h"
#include "MyConditionalPrior.h"

const double mas_to_rad = 4.84813681109536e-09;
extern const ComponentType component_type;
extern const bool use_per_antenna_jitters;
extern const bool use_single_jitter;
extern const bool use_offsets;

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
		
		// change the name of std::make_shared :)
		/**
		 * @brief Assign a prior distribution.
		 *
		 * This function defines, initializes, and assigns a prior distribution.
		 * Possible distributions are ...
		 *
		 * For example:
		 *
		 * @code{.cpp}
		 *          Cprior = make_prior<Uniform>(0, 1);
		 * @endcode
		 *
		 * @tparam T     ContinuousDistribution
		 * @tparam Args
		 * @param args   Arguments for constructor of distribution
		 * @return std::shared_ptr<T>
		*/
		template< class T, class... Args >
		std::shared_ptr<T> make_prior( Args&&... args ) { return std::make_shared<T>(args...); }

    private:
		bool fixed {false};
		int max_number_of_components {10};
		
		std::unordered_map<ComponentType, int> component_length{{circular, 4}, {sphere, 4}, {elliptical, 6}};
		
		std::vector<double> per_antenna_logjitter = std::vector<double>(Data::get_instance().n_antennas());
		std::shared_ptr<DNest4::ContinuousDistribution> Jsprior;
        double logjitter;
		std::shared_ptr<DNest4::ContinuousDistribution> Jprior;
		std::vector<double> per_antenna_offset = std::vector<double>(Data::get_instance().n_antennas());
		std::shared_ptr<DNest4::ContinuousDistribution> Oprior;
		unsigned int reference_antenna_number {4};
        unsigned int counter {0};
        DNest4::RJObject<MyConditionalPrior> components = DNest4::RJObject<MyConditionalPrior>(component_length[component_type], max_number_of_components, fixed, MyConditionalPrior());
        std::valarray<double> mu_real = std::valarray<double>(0., Data::get_instance().n_vis());
		std::valarray<double> mu_imag = std::valarray<double>(0., Data::get_instance().n_vis());
		// Dispersion for each point
		std::valarray<double> var = std::valarray<double>(0., Data::get_instance().n_vis());
		// Multiplicative amplitude offset for each point
		std::valarray<double> offset = std::valarray<double>(0., Data::get_instance().n_vis());
        void calculate_sky_mu();
		void calculate_var();
		void calculate_offset();
		void setPriors();
		
		// Pre-calculation of indexes
		std::vector<int> ant_ik;
		std::vector<int> ant_jk;
};

#endif //BAM_DNESTMODEL_H

#ifndef BAM_COMPONENT_H
#define BAM_COMPONENT_H

#include <cmath>
#include <valarray>
#include <vector>
#include <complex>
#include <iostream>
#include "RNG.h"
#include "DNest4.h"
#include "Utils.h"

const double mas_to_rad = 4.84813681109536e-09;


class Component {
    public:
		virtual ~Component() = default;
		virtual double get_logsize(double nu) = 0;
		virtual double get_logflux(double nu) = 0;
		virtual std::pair<double, double> get_pos(double nu) = 0;
		void ft(double nu, const std::valarray<double>& u, const std::valarray<double>& v);

        [[nodiscard]] std::valarray<double> get_mu_real() const
		{
            return mu_real;
        }

        [[nodiscard]] std::valarray<double> get_mu_imag() const
		{
            return mu_imag;
        }
		
        virtual void print(std::ostream& out) const = 0;
        virtual std::string description() const = 0;
        virtual void from_prior(DNest4::RNG& rng) = 0;
        virtual double perturb(DNest4::RNG& rng) = 0;
//        virtual double perturb_all_params(DNest4::RNG& rng) = 0;
//        virtual double perturb_position(DNest4::RNG& rng) = 0;
//        virtual double perturb_flux(DNest4::RNG& rng) = 0;
//        virtual double perturb_size(DNest4::RNG& rng) = 0;
//        virtual double perturb_flux_position(DNest4::RNG& rng) = 0;
//        virtual double perturb_size_position(DNest4::RNG& rng) = 0;
        // See also https://softwareengineering.stackexchange.com/a/337565 for unique_ptr
		// Should be implemented in each derived class:
        virtual Component* clone() = 0;
		enum Shape {Gaussian, Sphere};

    protected:
        // Contribution to the SkyModel prediction
        std::valarray<double> mu_real;
        std::valarray<double> mu_imag;
		Shape shape;
};


class JetComponent : public Component {
	public:
		JetComponent(Shape shape);
		JetComponent(const JetComponent& other);
		JetComponent* clone() override;
		double get_logsize(double nu) override;
		double get_logflux(double nu) override;
		std::pair<double, double> get_pos(double nu) override;
		void print(std::ostream& out) const override;
		std::string description() const override;
		void from_prior(DNest4::RNG& rng) override;
		double perturb(DNest4::RNG& rng) override;
	private:
		// RA, DEC
		double dx_{}, dy_{}, logsize_{}, lognu_max_{}, logS_max_{}, alpha_thick_{}, alpha_thin_{};
};


class CoreComponent : public Component {
	public:
		CoreComponent(Shape shape);
		CoreComponent(const CoreComponent& other);
		CoreComponent* clone() override;
		double get_logsize(double nu) override;
		double get_logflux(double nu) override;
		std::pair<double, double> get_pos(double nu) override;
		void print(std::ostream& out) const override;
		std::string description() const override;
		void from_prior(DNest4::RNG& rng) override;
		double perturb(DNest4::RNG& rng) override;
	private:
		// PA - from N to positive RA axis
//		double a_{}, PA_{}, logsize_1_{}, k_r_{}, logS_1_{}, alpha_{};
		double a_{}, PA_{}, logsize_1_{}, k_r_{}, lognu_max_{}, logS_max_{}, alpha_thick_{}, alpha_thin_{};
};

#endif //BAM_COMPONENT_H

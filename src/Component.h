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
		virtual double get_logsize(double nu) = 0;
		virtual double get_logflux(double nu) = 0;
		virtual std::pair<double, double> get_pos(double nu) = 0;
        virtual void ft(double nu, std::valarray<double> u, std::valarray<double> v) = 0;
//        virtual const size_t size() const = 0;

        [[nodiscard]] std::valarray<double> get_mu_real() const {
            return mu_real;
        }

        [[nodiscard]] std::valarray<double> get_mu_imag() const {
            return mu_imag;
        }

        virtual void print(std::ostream& out) const = 0;
        virtual std::string description() const = 0;
        virtual void from_prior(DNest4::RNG& rng) = 0;
        virtual double perturb(DNest4::RNG& rng) = 0;
        // See also https://softwareengineering.stackexchange.com/a/337565 for unique_ptr
		// Should be implemented in each derived class:
        virtual Component* clone() = 0;

    protected:
        // Contribution to the SkyModel prediction
        std::valarray<double> mu_real;
        std::valarray<double> mu_imag;
};


class GaussianComponent : public Component {
	public:
		GaussianComponent();
        void ft(double nu, std::valarray<double> u, std::valarray<double> v) override;
//		void print(std::ostream& out) const override {}
//        void from_prior(DNest4::RNG& rng) override {}
//        double perturb(DNest4::RNG& rng) override
//		{
//            std::cout << "Hope this is not called" << std::endl;
//            return 0.0;
//		}
};


class JetGaussianComponent : public GaussianComponent {
	public:
		JetGaussianComponent();
        JetGaussianComponent(const JetGaussianComponent& other);
        JetGaussianComponent* clone() override;
		double get_logsize(double nu) override;
		double get_logflux(double nu) override;
		std::pair<double, double> get_pos(double nu) override;
		void print(std::ostream& out) const override;
		std::string description() const override;
		void from_prior(DNest4::RNG& rng) override;
		double perturb(DNest4::RNG& rng) override;
	private:
		// RA, DEC
        double dx_, dy_, logsize_, lognu_max_, logS_max_, alpha_thick_, alpha_thin_;
};


class CoreGaussianComponent : public GaussianComponent {
	public:
		CoreGaussianComponent();
		CoreGaussianComponent(const CoreGaussianComponent& other);
		CoreGaussianComponent* clone() override;
		double get_logsize(double nu) override;
		double get_logflux(double nu) override;
		std::pair<double, double> get_pos(double nu) override;
		void print(std::ostream& out) const override;
		std::string description() const override;
		void from_prior(DNest4::RNG& rng) override;
		double perturb(DNest4::RNG& rng) override;
	private:
		// PA - from N to positive RA axis
		double a_, PA_, logsize_1_, k_r_, logS_1_, alpha_;
};

//class EGComponent : public Component {
//
//    public:
//        EGComponent();
//        void ft(double nu, std::valarray<double> u, std::valarray<double> v) override;
//        const size_t size() const override
//        {
//            return 9;
//        }
//        void print(std::ostream& out) const override;
//        void from_prior(DNest4::RNG& rng) override {}
//        double perturb(DNest4::RNG& rng) override {
//            std::cout << "Hope this is not called" << std::endl;
//            return 0.0; }
//
//    protected:
//        // Parameters of a single Gaussian
//        double dx_, dy_, logbmaj_, e_, bpa_, lognu_max, logS_max, alpha_thick, alpha_thin;
//};
//
//
//class CGComponent : public EGComponent {
//
//    public:
//        CGComponent();
//        void ft(double nu, std::valarray<double> u, std::valarray<double> v) override;
//        CGComponent(const CGComponent& other);
//        const size_t size() const override
//        {
//            return 7;
//        }
//        void print(std::ostream& out) const override;
//        std::string description() const override;
//        void from_prior(DNest4::RNG& rng) override;
//        double perturb(DNest4::RNG& rng) override;
//        CGComponent* clone() override;
//};
//
//
//class CoreCGComponent : public CGComponent {
//	public:
//		CoreCGComponent();
//		void ft(double nu, std::valarray<double> u, std::valarray<double> v) override;
//		CoreCGComponent(const CoreCGComponent& other);
//		const size_t size() const override
//		{
//			return 4;
//		}
//		void print(std::ostream& out) const override;
//		std::string description() const override;
//		void from_prior(DNest4::RNG& rng) override;
//		double perturb(DNest4::RNG& rng) override;
//		CoreCGComponent* clone() override;
//	private:
//		double a, PA, log_theta_1, k_r, logS_1, alpha;
//};


#endif //BAM_COMPONENT_H

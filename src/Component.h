#ifndef BAM_COMPONENT_H
#define BAM_COMPONENT_H

#include <cmath>
#include <valarray>
#include <vector>
#include <complex>
#include <iostream>
#include "RNG.h"
#include "Utils.h"

const double mas_to_rad = 4.84813681109536e-09;


class Component {
    public:
        virtual void ft(std::valarray<double> u, std::valarray<double> v) = 0;
        virtual const size_t size() const = 0;
        //virtual void set_param_vec(std::valarray<double> param) = 0;
        std::valarray<double> get_mu_real() const {
            return mu_real;
        }

        std::valarray<double> get_mu_imag() const {
            return mu_imag;
        }
        virtual void print(std::ostream& out) const = 0;
        virtual std::string description() const = 0;
        virtual void from_prior(DNest4::RNG& rng) = 0;
        virtual double perturb(DNest4::RNG& rng) = 0;
        // See also https://softwareengineering.stackexchange.com/a/337565 for unique_ptr
        virtual Component* clone() = 0;

    protected:
        // SkyModel prediction
        std::valarray<double> mu_real;
        std::valarray<double> mu_imag;
};


class DFComponent : public  Component {
    public:
        DFComponent();
        void ft(std::valarray<double> u, std::valarray<double> v) override;
        //void set_param_vec(std::valarray<double> param) override;
        const size_t size() const override
        {
            return 3;
        }
        void print(std::ostream& out) const override;
        void from_prior(DNest4::RNG& rng) override {}
        double perturb(DNest4::RNG& rng) override { return 0.0; }

    private:
        // Parameters of a single Delta Function
        double dx_, dy_, logflux_;
};


class EGComponent : public Component {

    public:
        EGComponent();
        void ft(std::valarray<double> u, std::valarray<double> v) override;
        //void set_param_vec(std::valarray<double> param) override;
        const size_t size() const override
        {
            return 6;
        }
        void print(std::ostream& out) const override;
        void from_prior(DNest4::RNG& rng) override {}
        double perturb(DNest4::RNG& rng) override {
            std::cout << "Hope this is not called" << std::endl;
            return 0.0; }

    protected:
        // Parameters of a single Gaussian
        double dx_, dy_, logflux_, logbmaj_, e_, bpa_;
};


class CGComponent : public EGComponent {

    public:
        CGComponent();
        CGComponent(const CGComponent& other);
        //void set_param_vec(std::valarray<double> param) override;
        const size_t size() const override
        {
            return 4;
        }
        void print(std::ostream& out) const override;
        std::string description() const override;
        void from_prior(DNest4::RNG& rng) override;
        double perturb(DNest4::RNG& rng) override;
        CGComponent* clone() override;
};


#endif //BAM_COMPONENT_H

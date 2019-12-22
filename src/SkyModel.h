#ifndef BSC_MODEL_H
#define BSC_MODEL_H

#include <vector>
#include "Component.h"
#include "RNG.h"


class SkyModel {
    public:
        SkyModel();
        SkyModel(const SkyModel& other);
        ~SkyModel();
        SkyModel&operator=(const SkyModel& other);

        void add_component(Component* component);
        void ft(const std::valarray<double>& u, const std::valarray<double>& v);
        void ft_from_all(const std::valarray<double>& u, const std::valarray<double>& v);
        std::valarray<double> get_mu_real() const { return mu_real; }
        std::valarray<double> get_mu_imag() const { return mu_imag; }
        void print(std::ostream& out) const;
        std::string description() const;
        void from_prior(DNest4::RNG& rng);
        // MH proposal for SkyModel. Returns logH
        double perturb(DNest4::RNG& rng);

    private:
        std::vector<Component*> components_;
        std::vector<bool> updated_;
        // SkyModel prediction
        std::valarray<double> mu_real;
        std::valarray<double> mu_imag;

};


#endif //BSC_MODEL_H

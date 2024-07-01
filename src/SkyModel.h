#ifndef BSC_MODEL_H
#define BSC_MODEL_H

#include <vector>
#include <unordered_map>
#include "Component.h"
#include "RNG.h"


class SkyModel
{
    public:
        SkyModel(size_t n_jet_components);
        SkyModel(const SkyModel& other);
        ~SkyModel();
        SkyModel&operator=(const SkyModel& other);
        SkyModel* clone();

        void add_component(Component* component);
		ArrayXcd ft_from_all(double nu, const ArrayXd& u, const ArrayXd& v);
        void print(std::ostream& out) const;
		std::string print() const;
        [[nodiscard]] std::string description() const;
        void from_prior(DNest4::RNG& rng);
        // MH proposal for SkyModel. Returns logH
        double perturb(DNest4::RNG& rng);
		std::pair<double, double> get_core_position(double nu);

    private:
        std::vector<Component*> components_{};
};


#endif //BSC_MODEL_H

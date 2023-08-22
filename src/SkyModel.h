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
		ArrayXcd ft_from_perturbed(double nu, const ArrayXd& u, const ArrayXd& v);
//        ArrayXcd get_mu() const { return mu; }
		void set_perturbed(std::vector<bool> perturbed);
		std::vector<bool> get_perturbed();
		void reset_perturbed();
        void print(std::ostream& out) const;
        [[nodiscard]] std::string description() const;
        void from_prior(DNest4::RNG& rng);
        // MH proposal for SkyModel. Returns logH
        double perturb(DNest4::RNG& rng);

    private:
        std::vector<bool> perturbed{};
        std::vector<Component*> components_{};
        // SkyModel prediction
//		ArrayXcd mu{};
};


#endif //BSC_MODEL_H

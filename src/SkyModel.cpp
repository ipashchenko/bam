#include <SkyModel.h>
#include <iostream>
#include "Data.h"

SkyModel::SkyModel() = default;


SkyModel::SkyModel(const SkyModel &other) {
    for (auto other_comp : other.components_) {
        Component* comp = other_comp->clone();
        components_.push_back(comp);
    }
    mu_real = other.mu_real;
    mu_imag = other.mu_imag;
}


SkyModel& SkyModel::operator=(const SkyModel& other) {
    if (&other == this) {
        return *this;
    }

    for (auto comp : components_) {
        delete comp;
    }
    components_.clear();

    for (auto other_comp : other.components_) {
        Component* comp = other_comp->clone();
        components_.push_back(comp);
    }
    mu_real = other.mu_real;
    mu_imag = other.mu_imag;
    return *this;
}



SkyModel::~SkyModel() {
    for (auto comp : components_) {
        delete comp;
    }
    components_.clear();
}


void SkyModel::add_component(Component *component) {
    components_.push_back(component);
}

void SkyModel::ft(const std::valarray<double>& u, const std::valarray<double>& v) {
    // Traverse components and sum differences of new and old predictions for updated components.
    for (auto comp : components_) {
        if(comp->is_updated) {
            comp->ft(u, v);
            // TODO: Here I can add to ``mu_real``, ``mu_imag`` w/o declaring new arrays. However ``mu_real/imag`` must
            // be already initialized. They are initialized to zeros in ``SkyModel.from_prior`` and first
            // ``ft_from_all`` fills them with predictions.
            mu_real += (comp->get_mu_real() - comp->get_mu_real_old());
            mu_imag += (comp->get_mu_imag() - comp->get_mu_imag_old());
            comp->is_updated = false;
            comp->update_old();
            break;
        }
    }
}


void SkyModel::ft_from_all(const std::valarray<double>& u, const std::valarray<double>& v) {
    // Zero prediction
    mu_real *= 0.0;
    mu_imag *= 0.0;

    for (auto comp : components_) {
        comp->ft(u, v);
        // TODO: Here I can add to ``mu_real``, ``mu_imag`` w/o declaring new arrays
        mu_real += comp->get_mu_real();
        mu_imag += comp->get_mu_imag();
        if(comp->is_updated) {
            comp->is_updated = false;
            comp->update_old();
        }
    }
}


void SkyModel::print(std::ostream &out) const
{
    for (auto comp: components_) {
        comp->print(out);
    }
}


void SkyModel::from_prior(DNest4::RNG &rng) {
    const std::valarray<double>& u = Data::get_instance().get_u();
    std::valarray<double> zero (0.0, u.size());
    mu_real = zero;
    mu_imag = zero;
    for (auto comp: components_) {
        comp->from_prior(rng);
        comp->is_updated = true;
    }
}


double SkyModel::perturb(DNest4::RNG &rng) {
    int which = rng.rand_int(components_.size());
    components_[which]->is_updated = true;
    return components_[which]->perturb(rng);
}


std::string SkyModel::description() const {
    std::string descr;
    for (auto comp: components_) {
        descr += comp->description();
        descr += " ";
    }
    descr.pop_back();
    return descr;
}
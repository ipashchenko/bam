#include <SkyModel.h>
#include <iostream>


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
    std::valarray<double> real (0.0, u.size());
    std::valarray<double> imag (0.0, u.size());
    for (auto comp : components_) {
        comp->ft(u, v);
        real = real + comp->get_mu_real();
        imag = imag + comp->get_mu_imag();
    }
    mu_real = real;
    mu_imag = imag;
}


void SkyModel::print(std::ostream &out) const
{
    for (auto comp: components_) {
        comp->print(out);
    }
}


void SkyModel::from_prior(DNest4::RNG &rng) {
    for (auto comp: components_) {
        comp->from_prior(rng);
    }
}


double SkyModel::perturb(DNest4::RNG &rng) {
    int which = rng.rand_int(components_.size());
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
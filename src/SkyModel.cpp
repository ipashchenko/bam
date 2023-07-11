#include <SkyModel.h>
#include <iostream>
#include "Data.h"

SkyModel::SkyModel() = default;


SkyModel::SkyModel(const SkyModel &other)
{
    for (auto other_comp : other.components_)
	{
        Component* comp = other_comp->clone();
        components_.push_back(comp);
    }
    mu_real = other.mu_real;
    mu_imag = other.mu_imag;
}


SkyModel& SkyModel::operator=(const SkyModel& other)
{
    if (&other == this)
	{
        return *this;
    }

    for (auto comp : components_)
	{
        delete comp;
    }
    components_.clear();

    for (auto other_comp : other.components_)
	{
        Component* comp = other_comp->clone();
        components_.push_back(comp);
    }
    mu_real = other.mu_real;
    mu_imag = other.mu_imag;
    return *this;
}



SkyModel::~SkyModel()
{
    for (auto comp : components_)
	{
        delete comp;
    }
    components_.clear();
}


void SkyModel::add_component(Component *component)
{
    components_.push_back(component);
}

void SkyModel::ft_from_all(double nu, const std::valarray<double>& u, const std::valarray<double>& v)
{
    // Zero prediction
	std::valarray<double> zero (0.0, u.size());
    mu_real = zero;
    mu_imag = zero;

    for (auto comp : components_)
	{
        comp->ft(nu, u, v);
        mu_real += comp->get_mu_real();
        mu_imag += comp->get_mu_imag();
    }
}


void SkyModel::print(std::ostream &out) const
{
    for (auto comp: components_)
	{
        comp->print(out);
    }
}


void SkyModel::from_prior(DNest4::RNG &rng)
{
    for (auto comp: components_)
	{
        comp->from_prior(rng);
    }
}


double SkyModel::perturb(DNest4::RNG &rng) {
    int which = rng.rand_int(components_.size());
//    components_[which]->is_updated = true;
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
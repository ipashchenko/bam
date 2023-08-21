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
    for (auto other_comp : other.old_components_)
	{
        Component* comp = other_comp->clone();
        old_components_.push_back(comp);
    }
    updated_ = other.updated_;
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

    for (auto comp : old_components_)
	{
        delete comp;
    }
    old_components_.clear();

    for (auto other_comp : other.components_)
	{
        Component* comp = other_comp->clone();
        components_.push_back(comp);
    }
    for (auto other_comp : other.old_components_)
	{
        Component* comp = other_comp->clone();
        old_components_.push_back(comp);
    }
    updated_ = other.updated_;
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
    for (auto comp : old_components_)
	{
        delete comp;
    }
    old_components_.clear();
}


void SkyModel::add_component(Component *component)
{
    components_.push_back(component);
    updated_.push_back(true);
    Component* old_component = component->clone();
    old_components_.push_back(old_component);
}

void SkyModel::ft_from_all(double nu, const std::valarray<double>& u, const std::valarray<double>& v)
{
    // Zero prediction
	std::valarray<double> zero (0.0, u.size());
    mu_real = zero;
    mu_imag = zero;

    for (size_t i=0; i < components_.size(); i++)
	{   auto comp = components_[i];
        comp->ft(nu, u, v);
        mu_real += comp->get_mu_real();
        mu_imag += comp->get_mu_imag();
		updated_[i] = false;
    }
}


void SkyModel::ft(double nu, const std::valarray<double>& u, const std::valarray<double>& v)
{
	for (size_t i=0; i < components_.size(); i++)
	{
        auto comp = components_[i];
		auto old_comp = old_components_[i];
        if(updated_[i])
		{
			comp->ft(nu, u, v);
			mu_real += (comp->get_mu_real() - old_comp->get_mu_real());
			mu_imag += (comp->get_mu_imag() - old_comp->get_mu_imag());
			updated_[i] = false;
		}
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
    for (size_t i=0; i < components_.size(); i++)
	{
        auto comp = components_[i];
        comp->from_prior(rng);
        updated_[i] = true;
        Component* old_comp = comp->clone();
        old_components_[i] = old_comp;
    }
}


double SkyModel::perturb(DNest4::RNG &rng) {
    int which = rng.rand_int(components_.size());
    updated_[which] = true;
    return components_[which]->perturb(rng);
}


std::string SkyModel::description() const {
    std::string descr;
    for (auto comp: components_) {
        descr += comp->description();
        descr += "\t";
    }
    descr.pop_back();
    return descr;
}
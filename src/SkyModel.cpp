#include <SkyModel.h>
#include <iostream>
#include "Data.h"


SkyModel::SkyModel(size_t n_jet_components) 
{
	auto* comp = new CoreComponent(Component::Gaussian);
	this->add_component(comp);
    for (int i=0; i<n_jet_components; i++) {
        auto* comp = new JetComponent(Component::Gaussian);
        this->add_component(comp);
    }
}


SkyModel::SkyModel(const SkyModel &other)
{
    for (auto other_comp : other.components_)
	{
        Component* comp = other_comp->clone();
        components_.push_back(comp);
    }
    perturbed = other.perturbed;
    mu = other.mu;
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
    perturbed = other.perturbed;
    mu = other.mu;
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


SkyModel* SkyModel::clone()
{
    return new SkyModel(*this);
}


void SkyModel::set_perturbed(std::vector<bool> new_perturbed)
{
	perturbed = new_perturbed;
}

std::vector<bool> SkyModel::get_perturbed()
{
	return perturbed;
}

void SkyModel::add_component(Component *component)
{
    components_.push_back(component);
    perturbed.push_back(false);
}

void SkyModel::ft_from_all(double nu, const ArrayXd& u, const ArrayXd& v)
{
    // Zero prediction
	ArrayXcd zero = ArrayXcd::Zero(u.size());
    mu = zero;

    for (auto comp : components_)
	{
		mu += comp->ft(nu, u, v);
    }
}


void SkyModel::ft_from_perturbed(double nu, const ArrayXd& u, const ArrayXd& v)
{
	// Zero prediction
	mu *= 0.0;
    for(size_t i=0; i<perturbed.size(); i++)
    {
        if(perturbed[i])
        {
            auto comp = components_[i];
            mu += comp->ft(nu, u, v);
            perturbed[i] = false;
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
    for (auto comp : components_)
	{
        comp->from_prior(rng);
    }
}


double SkyModel::perturb(DNest4::RNG &rng) {
    int which = rng.rand_int(components_.size());
    perturbed[which] = true;
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
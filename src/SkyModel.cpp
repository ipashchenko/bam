#include <SkyModel.h>
#include <iostream>
#include "Data.h"

#ifdef NDEBUG
#define DEBUG(x)
#else
#define DEBUG(x) do { std::cerr << x << std::endl; } while (0)
#endif


SkyModel::SkyModel(size_t n_jet_components) 
{
	auto* comp = new CoreComponent(Component::Gaussian);
	DEBUG("Added core component");
	this->add_component(comp);
    for (int i=0; i<n_jet_components; i++) {
        auto* comp = new JetComponent(Component::Gaussian);
        this->add_component(comp);
		DEBUG("Added jet component : " + comp->print());
    }
}


SkyModel::SkyModel(const SkyModel &other)
{
    for (auto other_comp : other.components_)
	{
        Component* comp = other_comp->clone();
        components_.push_back(comp);
    }
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


void SkyModel::add_component(Component *component)
{
    components_.push_back(component);
}

ArrayXcd SkyModel::ft_from_all(double nu, const ArrayXd& u, const ArrayXd& v)
{
    // Zero prediction
    ArrayXcd mu = ArrayXcd::Zero(u.size());

    for (auto comp : components_)
	{
		DEBUG("SM.ft_from_all : mu[0] = " + std::to_string(mu[0].real()) + ", " + std::to_string(mu[0].imag()));
		mu = (mu + comp->ft(nu, u, v)).eval();
    }
	
	return mu;
}

void SkyModel::print(std::ostream &out) const
{
    for (auto comp: components_)
	{
        comp->print(out);
    }
}

std::string SkyModel::print() const
{
	std::string out("SkyModel print : \n");
	for (auto comp: components_)
	{
		out += comp->print();
	}
	return out;
}

void SkyModel::from_prior(DNest4::RNG &rng)
{
    for (auto comp : components_)
	{
        comp->from_prior(rng);
    }
}


double SkyModel::perturb(DNest4::RNG &rng) {
    int which = rng.rand_int(2*components_.size());
	DEBUG("which = " + std::to_string(which));
	// More often update core component
	if(which >= components_.size())
	{
		which = 0;
	}
	DEBUG("In SkyModel.perturb() - perturbing component # " + std::to_string(which));
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


std::pair<double, double> SkyModel::get_core_position(double nu)
{
	return components_[0]->get_pos(nu);
}

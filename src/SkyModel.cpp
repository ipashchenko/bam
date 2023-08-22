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
		DEBUG("Added jet component");
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

void SkyModel::reset_perturbed()
{
	fill(perturbed.begin(), perturbed.end(), false);
}

void SkyModel::add_component(Component *component)
{
    components_.push_back(component);
    perturbed.push_back(false);
}

void SkyModel::ft_from_all(double nu, const ArrayXd& u, const ArrayXd& v)
{
    // Zero prediction
    mu = ArrayXcd::Zero(u.size());

    for (auto comp : components_)
	{
		mu += comp->ft(nu, u, v);
    }
}


void SkyModel::ft_from_perturbed(double nu, const ArrayXd& u, const ArrayXd& v)
{
	if (std::find(std::begin(perturbed), std::end(perturbed), true) == std::end(perturbed)) // All false
	{
		std::cout << "All perturbed = false, but we in SkyModel.ft_from_perturbed!!!!!!\n";
	}
	// Zero prediction
	// TODO: Can't I multiply mu on ``0``?
	mu = ArrayXcd::Zero(u.size());
    for(size_t i=0; i<perturbed.size(); i++)
    {
        if(perturbed[i])
        {
			DEBUG("FT only from component # " + std::to_string(i));
            auto comp = components_[i];
            mu = (mu + comp->ft(nu, u, v)).eval();
			// FIXME: I can't reset perturbed flag here, because it calculates FT for several bands!
//            perturbed[i] = false;
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
	DEBUG("In SkyModel.perturb() - perturbing component # " + std::to_string(which));
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
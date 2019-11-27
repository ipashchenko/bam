#include <iostream>
#include <Distributions/ContinuousDistribution.h>
#include "DNestModel.h"
#include "Data.h"
#include "Component.h"


DNestModel::DNestModel() {
    sky_model = new SkyModel();
    int ncomp = 6;
    for (int i=0; i<ncomp; i++) {
        auto* comp = new CGComponent();
        sky_model->add_component(comp);
    }
}


DNestModel::~DNestModel() {
    delete sky_model;
}


DNestModel::DNestModel(const DNestModel& other) {
    sky_model = new SkyModel(*other.sky_model);
    mu_real = other.mu_real;
    mu_imag = other.mu_imag;
}


DNestModel& DNestModel::operator=(const DNestModel& other) {
    if (this != &other) {
        *(sky_model) = *(other.sky_model);
        mu_real = other.mu_real;
        mu_imag = other.mu_imag;
    }
    return *this;
}


void DNestModel::from_prior(DNest4::RNG &rng) {
    // I do this because in ``calculate_sky_mu`` ``mu_real`` and ``mu_imag`` are multiplied and added.
    const std::valarray<double>& u = Data::get_instance().get_u();
    std::valarray<double> zero (0.0, u.size());
    mu_real = zero;
    mu_imag = zero;

    sky_model->from_prior(rng);
    calculate_sky_mu(false);
}


double DNestModel::perturb(DNest4::RNG &rng) {
    double logH = 0.;

    logH += sky_model->perturb(rng);

    // Pre-reject
    if(rng.rand() >= exp(logH)) {
        return -1E300;
    }
    else
        logH = 0.0;

    // This shouldn't be called in case of pre-rejection
    calculate_sky_mu(false);
    return logH;
}


void DNestModel::calculate_sky_mu(bool update) {
    const std::valarray<double>& u = Data::get_instance().get_u();
    const std::valarray<double>& v = Data::get_instance().get_v();

    // FT (calculate SkyModel prediction)
    sky_model->ft(u, v);
    mu_real = sky_model->get_mu_real();
    mu_imag = sky_model->get_mu_imag();
}


double DNestModel::log_likelihood() const {
    const std::valarray<double>& vis_real = Data::get_instance().get_vis_real();
    const std::valarray<double>& vis_imag = Data::get_instance().get_vis_imag();
    const std::valarray<double>& sigma = Data::get_instance().get_sigma();

    // Variance
    const std::valarray<double> var = sigma*sigma;
    // Complex Gaussian sampling distribution
    std::valarray<double> result = -log(2*M_PI*var) - 0.5*(pow(vis_real - mu_real, 2) +
        pow(vis_imag - mu_imag, 2))/var;
    double loglik = result.sum();
    return loglik;

}


void DNestModel::print(std::ostream &out) const {
    sky_model->print(out);
}


std::string DNestModel::description() const
{
    std::string descr;
    descr += sky_model->description();
    return descr;
}

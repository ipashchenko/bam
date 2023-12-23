#include <iostream>
#include <Distributions/ContinuousDistribution.h>
#include "DNestModel.h"
#include "Data.h"

// TODO: Uniform priors on position to bracket!
DNestModel::DNestModel() :
    logjitter(0.0),
    counter(0),
    components(4, 20, false, MyConditionalPrior(-2., 20., -2., 2.), DNest4::PriorType::log_uniform) {}


void DNestModel::from_prior(DNest4::RNG &rng) {
    // I do this because in ``calculate_sky_mu`` ``mu_real`` and ``mu_imag`` are multiplied and added.
    const std::valarray<double>& u = Data::get_instance().get_u();
    std::valarray<double> zero (0.0, u.size());
    mu_real = zero;
    mu_imag = zero;
    logjitter = -4.0 + 2.0*rng.randn();
    components.from_prior(rng);
    components.consolidate_diff();
    calculate_sky_mu();
}


double DNestModel::perturb(DNest4::RNG &rng) {
    double logH = 0.;

    // Perturb jitter
    if(rng.rand() <= 0.1) {
        logH -= -0.5*pow((logjitter+4)/2.0, 2.0);
        logjitter += 2.0*rng.randh();
        logH += -0.5*pow((logjitter+4)/2.0, 2.0);

        // Pre-reject
        if(rng.rand() >= exp(logH)) {
            return -1E300;
        }
        else
            logH = 0.0;
        // No need to re-calculate model. Just calculate loglike.
    }

    // Perturb SkyModel
    else {
        logH += components.perturb(rng);
        components.consolidate_diff();
        // After this ``removed`` is empty and gone to ``added`` with ``-`` sign. We use ``added`` when
        // ``update = True``. Else we use all components.

        // Pre-reject
        if(rng.rand() >= exp(logH)) {
            return -1E300;
        }
        else
            logH = 0.0;

        // This shouldn't be called in case of pre-rejection
        calculate_sky_mu();
    }
    return logH;
}


void DNestModel::calculate_sky_mu() {
    const std::valarray<double>& u = Data::get_instance().get_u();
    const std::valarray<double>& v = Data::get_instance().get_v();

    std::valarray<double> theta;
    double c;
    std::valarray<double> b;
    std::valarray<double> ft;

    // Update or from scratch?
    bool update = (components.get_added().size() < components.get_components().size()) &&
        (counter <= 20);
    // Get the components
    const std::vector<std::vector<double>>& comps = (update)?(components.get_added()):
                                                    (components.get_components());

    if(!update)
    {
        // Zero the sky model prediction
        mu_real *= 0.0;
        mu_imag *= 0.0;
        //std::cout << "full " << counter << "\n";
        counter = 0;

    } else {
        //std::cout << "diff " << counter << "\n";
        counter++;
    }

    for(const auto& comp: comps)
    {
        // Phase shift due to not being in a phase center
        theta = 2*M_PI*mas_to_rad*(u*comp[0]+v*comp[1]);
        // Calculate FT of a Gaussian in a phase center
        c = pow(M_PI*exp(comp[3])*mas_to_rad, 2)/(4.*log(2.));
        ft = comp[2] * exp(-c*(u*u + v*v));
        // Prediction of visibilities
        mu_real += ft*cos(theta);
        mu_imag += ft*sin(theta);
    }
}


double DNestModel::log_likelihood() const {
    const std::valarray<double>& vis_real = Data::get_instance().get_vis_real();
    const std::valarray<double>& vis_imag = Data::get_instance().get_vis_imag();
    const std::valarray<double>& sigma = Data::get_instance().get_sigma();

    // Variance
    // TODO: Keep variance in data not to square sigma each time
    const std::valarray<double> var = sigma*sigma;
    // Complex Gaussian sampling distribution
    std::valarray<double> result = -log(2*M_PI*(var+exp(2.0*logjitter))) - 0.5*(pow(vis_real - mu_real, 2) +
        pow(vis_imag - mu_imag, 2))/(var+exp(2.0*logjitter))   ;
    double loglik = result.sum();
    return loglik;

}


void DNestModel::print(std::ostream &out) const {
    out << logjitter << '\t';
    components.print(out);
}


std::string DNestModel::description() const
{
    std::string descr;

    // Anything printed by DNestModel::print (except the last line)
    descr += "logjitter ";

    // The rest is all what happens when you call .print on an RJObject
    descr += " dim_components max_num_components ";

    // Then the hyperparameters (i.e. whatever MyConditionalPrior::print prints)
    descr += " typical_flux dev_log_flux typical_radius dev_log_radius ";

    // Then the actual number of components
    descr += " num_components ";

    // Then it's all the components, padded with zeros
    // max_num_components is known in this model, so that's how far the
    // zero padding goes.
    for(int i=0; i<components.get_max_num_components(); ++i)
        descr += " x[" + std::to_string(i) + "] ";
    for(int i=0; i<components.get_max_num_components(); ++i)
        descr += " y[" + std::to_string(i) + "] ";
    for(int i=0; i<components.get_max_num_components(); ++i)
        descr += " logflux[" + std::to_string(i) + "] ";
    for(int i=0; i<components.get_max_num_components(); ++i)
        descr += " logbmaj[" + std::to_string(i) + "] ";

    return descr;
}

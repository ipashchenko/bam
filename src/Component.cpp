#include "Component.h"
#include "Data.h"
#include <cmath>
#include <iostream>
#include <RNG.h>


DFComponent::DFComponent() : dx_(0.0), dy_(0.0), logflux_(0.0)
{

}


void DFComponent::ft(std::valarray<double> u, std::valarray<double> v)
{
    std::valarray<double> theta;
    // Phase shift due to not being in a phase center
    theta = 2*M_PI*mas_to_rad*(u*dx_+v*dy_);

    // Prediction of visibilities
    mu_real = exp(logflux_)*cos(theta);
    mu_imag = exp(logflux_)*sin(theta);
}


void DFComponent::print(std::ostream &out) const
{
    out << "x=" << dx_ << ", " << "y=" << dy_ << ", " << "logflux=" << logflux_ << "\n";

}


EGComponent::EGComponent() : dx_(0.0), dy_(0.0), logflux_(0.0), logbmaj_(0.0), e_(1.0), bpa_(0.0)
{

}


void EGComponent::ft(std::valarray<double> u, std::valarray<double> v)
{
    std::valarray<double> theta;
    double c;
    std::valarray<double> b;
    std::valarray<double> ft;

    // Phase shift due to not being in a phase center
    theta = 2*M_PI*mas_to_rad*(u*dx_+v*dy_);
    // Calculate FT of a Gaussian in a phase center
    c = pow(M_PI*exp(logbmaj_)*mas_to_rad, 2)/(4.*log(2.));
    b = e_*e_*pow((u*cos(bpa_) - v*sin(bpa_)), 2) + pow((u*sin(bpa_)+v*cos(bpa_)), 2);
    ft = exp(logflux_-c*b);

    // Prediction of visibilities
    mu_real = ft*cos(theta);
    mu_imag = ft*sin(theta);
}


void EGComponent::print(std::ostream &out) const
{
    out << "x=" << dx_ << ", " << "y=" << dy_ << ", " << "logflux=" << logflux_ << ", "<< "logsize=" << logbmaj_ << ", "<< "e=" << e_ << ", "<< "bpa=" << bpa_ << "\n";

}


CGComponent::CGComponent() : EGComponent() {}


void CGComponent::ft(std::valarray<double> u, std::valarray<double> v)
{
    std::valarray<double> theta;
    double c;
    //std::valarray<double> b;
    std::valarray<double> ft;

    // Phase shift due to not being in a phase center
    theta = 2*M_PI*mas_to_rad*(u*dx_+v*dy_);
    // Calculate FT of a Gaussian in a phase center
    c = pow(M_PI*exp(logbmaj_)*mas_to_rad, 2)/(4.*log(2.));
    //b = e_*e_*pow((u*cos(bpa_) - v*sin(bpa_)), 2) + pow((u*sin(bpa_)+v*cos(bpa_)), 2);
    // TODO: Keep u*u and v*v already computed
    ft = exp(logflux_-c*(u*u + v*v));

    // Prediction of visibilities
    mu_real = ft*cos(theta);
    mu_imag = ft*sin(theta);
}


void CGComponent::print(std::ostream &out) const
{
    out << dx_ << '\t' << dy_ << '\t' << logflux_ << '\t' << logbmaj_ << '\t';
}


void CGComponent::from_prior(DNest4::RNG &rng) {
     const std::valarray<double>& u = Data::get_instance().get_u();
     std::valarray<double> zero (0.0, u.size());
     mu_real = zero;
     mu_imag = zero;
     mu_real_old = zero;
     mu_imag_old = zero;
     // Normal diffuse prior for x & y
     dx_ = 5.0*rng.randn();
     dy_ = 5.0*rng.randn();
     // Log-normal prior for flux and bmaj
     logflux_ = -2.0 + 1.0*rng.randn();
     logbmaj_ = -1.0 + 2.0*rng.randn();
    // Dependent priors - from MOJAVE 15 GHz modelfits.
    //logflux_ = -1.74 + 1.62*rng.randn();
    //logbmaj_ = -0.71 - 0.44*logflux_ + 0.84*rng.randn();
    // Dependent priors - modified for 43 GHz.
//    logflux_ = -2.29 + 1.62*rng.randn();
//    logbmaj_ = -1.71 - 0.44*logflux_ + 0.84*rng.randn();
}


double CGComponent::perturb(DNest4::RNG &rng) {
    double log_H = 0.;
    int which = rng.rand_int(4);
    if(which%4 == 0)
    {
        log_H -= -0.5*pow(dx_/5.0, 2);
        dx_ += 5.0*rng.randh();
        log_H += -0.5*pow(dx_/5.0, 2);
    }
    else if(which%4 == 1)
    {
        log_H -= -0.5*pow(dy_/5.0, 2);
        dy_ += 5.0*rng.randh();
        log_H += -0.5*pow(dy_/5.0, 2);
    }
    else if(which%4 == 2)
    {
        log_H -= -0.5*pow((logflux_+2.0)/1.0, 2);
        logflux_ += 1.0*rng.randh();
        log_H += -0.5*pow((logflux_+2.0)/1.0, 2);

//        log_H -= -0.5*pow((logflux_+1.74)/1.62, 2);
//        logflux_ += 1.62*rng.randh();
//        log_H += -0.5*pow((logflux_+1.74)/1.62, 2);
    }
    else
    {
        log_H -= -0.5*pow((logbmaj_+1.0)/2.0, 2);
        logbmaj_ += 2.0*rng.randh();
        log_H += -0.5*pow((logbmaj_+1.0)/2.0, 2);

//        log_H -= -0.5*pow((logbmaj_+0.71+0.44*logflux_)/0.84, 2);
//        logbmaj_ += 0.84*rng.randh();
//        log_H += -0.5*pow((logbmaj_+0.71+0.44*logflux_)/0.84, 2);
    }
    return log_H;
}


CGComponent *CGComponent::clone() {
    return new CGComponent(*this);
}


CGComponent::CGComponent(const CGComponent &other) {
    dx_ = other.dx_;
    dy_ = other.dy_;
    logflux_ = other.logflux_;
    logbmaj_ = other.logbmaj_;
    mu_real = other.get_mu_real();
    mu_imag = other.get_mu_imag();
    mu_real_old = other.get_mu_real_old();
    mu_imag_old = other.get_mu_imag_old();
}


std::string CGComponent::description() const {
    std::string descr;
    descr += "x y flux bmaj";
    return descr;
}


#ifndef BSC_SRC_LOGNORMAL_H_
#define BSC_SRC_LOGNORMAL_H_

#include <Distributions/Gaussian.h>
#include "Distributions/ContinuousDistribution.h"


class LogNormal : public DNest4::ContinuousDistribution {
    private:
        DNest4::Gaussian gauss;
    public:
        LogNormal(double center=0.0, double width=1.0);
        double generate(DNest4::RNG& rng) const ;
        double cdf(double x) const;
        double cdf_inverse(double x) const;
        double log_pdf(double x) const;
};

#endif //BSC_SRC_LOGNORMAL_H_

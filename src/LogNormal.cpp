#include "Distributions/Gaussian.h"
#include "LogNormal.h"


LogNormal::LogNormal(double center, double width) :
    gauss(center, width)
{ }


double LogNormal::cdf(double x) const
{
    return gauss.cdf(log(x));
}


double LogNormal::cdf_inverse(double x) const
{
    return exp(gauss.cdf_inverse(x));
}


double LogNormal::log_pdf(double x) const
{
    return gauss.log_pdf(log(x));
}


double LogNormal::generate(DNest4::RNG &rng) const {
    return exp(gauss.generate(rng));
}

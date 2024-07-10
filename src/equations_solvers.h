#ifndef MHD_SOLVER_INCLUDE_EQUATIONS_SOLVERS_H
#define MHD_SOLVER_INCLUDE_EQUATIONS_SOLVERS_H

#include <iostream>
#include <tuple>
#include <boost/math/tools/roots.hpp>
#include <limits>
#include "utils.h"


// This aimed to solve n_local given M_local, n_ext, c_ext and eta
template <class T>
struct functor_n1deriv
{
        // Functor returning both 1st derivative.
        functor_n1deriv(T const& Msq, T const& n_ext, T const& cs_to_c_sq_ext, T const& eta, T const& Gamma) :
                Msq_(Msq), n_ext_(n_ext), cs_to_c_sq_ext_(cs_to_c_sq_ext), eta_(eta), Gamma_(Gamma)
        {}
        std::tuple<T, T> operator()(T const& x)
        {
            // Return both f(x) and f'(x).
            T fx = Msq_ / (4 * M_PI * eta_ * eta_ * m_e * c * c) * x - 1./(Gamma_ - 1.0)*cs_to_c_sq_ext_*pow(x/n_ext_, Gamma_ - 1.0) - 1.0;
            T dx = Msq_ / (4 * M_PI * eta_ * eta_ * m_e * c * c) - 1. / (pow(n_ext_, Gamma_ - 1.0)) * cs_to_c_sq_ext_ * pow(x, Gamma_ - 2.0);
            return std::make_tuple(fx, dx);  // 'return' fx, dx.
        }
    private:
        T Gamma_;
        T Msq_;
        T n_ext_;
        T cs_to_c_sq_ext_;
        T eta_;
};


template <class T>
struct functor_1deriv
{
        // Functor returning both 1st derivative.
        functor_1deriv(T const& Gamma, T const& cs_sq_0, T const& Msq_ratio) : Gamma_(Gamma), cs_sq_0_(cs_sq_0), Msq_ratio_(Msq_ratio)
        {}
        std::tuple<T, T> operator()(T const& x)
        {
            // Return both f(x) and f'(x).
            T fx = (1.0+cs_sq_0_/(Gamma_-1.0))*(Msq_ratio_)*x - cs_sq_0_*pow(x, Gamma_-1.0)/(Gamma_-1.0) - 1.0;
            T dx = (1.0+cs_sq_0_/(Gamma_-1.0))*(Msq_ratio_) - cs_sq_0_*pow(x, Gamma_-2.0);
            return std::make_tuple(fx, dx);  // 'return' fx, dx.
        }
    private:
        T Gamma_;
        // Squared ratio of sound speed (at the jet boundary) to speed of light
        T cs_sq_0_;
        // Squared ratio of the local Mach speed to Mach speed at the jet boundary
        T Msq_ratio_;
};


template <class T>
struct functor_2deriv
{
        // Functor returning both 1st and 2nd derivatives.
        functor_2deriv(T const& Gamma, T const& cs_sq_0, T const& Msq_ratio) : Gamma_(Gamma), cs_sq_0_(cs_sq_0), Msq_ratio_(Msq_ratio)
        {}
        std::tuple<T, T, T> operator()(T const& x)
        {
            // Return both f(x) and f'(x) and f''(x).
            T fx = (1.0+cs_sq_0_/(Gamma_-1.0))*(Msq_ratio_)*x - cs_sq_0_*pow(x, Gamma_-1.0)/(Gamma_-1.0) - 1.0;
            T dx = (1.0+cs_sq_0_/(Gamma_-1.0))*(Msq_ratio_) - cs_sq_0_*pow(x, Gamma_-2.0);
            T d2x = -cs_sq_0_*(Gamma_-2.0)*pow(x, Gamma_-3.0);;
            return std::make_tuple(fx, dx, d2x);  // 'return' fx, dx and d2x.
        }
    private:
        T Gamma_;
        // Squared ratio of sound speed (at the jet boundary) to speed of light
        T cs_sq_0_;
        // Squared ratio of the local Mach speed to Mach speed at the jet boundary
        T Msq_ratio_;
};

template <class T>
T find_n_to_n0_2deriv(T Gamma, T cs_sq_0, T Msq_ratio)
{
    using namespace boost::math::tools;
    const int digits = std::numeric_limits<T>::digits;  // Maximum possible binary digits accuracy for type T.
    // digits used to control how accurate to try to make the result.
    int get_digits = static_cast<int>(digits * 0.4);    // Accuracy triples with each step, so stop when just
    //std::cout << "n digits = " << get_digits << "\n";
    // over one third of the digits are correct.
    boost::uintmax_t maxit = 30;
    //std::cout << "Gamma = " << Gamma << ", cs_sq_0 = " << cs_sq_0 << ", Msq_ratio = " << Msq_ratio << "\n";
    //T result = halley_iterate(functor_2deriv<T>(Gamma, cs_sq_0, Msq_ratio), 0.1, 0.0, 1E+05, get_digits, maxit);
    T result = newton_raphson_iterate(functor_1deriv<T>(Gamma, cs_sq_0, Msq_ratio), 0.1, 0.0, 1E+05, get_digits, maxit);
    //std::cout << "Found res = " << result << "\n";
    return result;
}


template <class T>
T find_n_local_1deriv(T Msq_local, T n_ext, T cs_to_c_sq_ext, T eta, T Gamma)
{
    using namespace boost::math::tools;
    const int digits = std::numeric_limits<T>::digits;  // Maximum possible binary digits accuracy for type T.
    // digits used to control how accurate to try to make the result.
    int get_digits = static_cast<int>(digits * 0.4);    // Accuracy triples with each step, so stop when just
    //std::cout << "n digits = " << get_digits << "\n";
    // over one third of the digits are correct.
    boost::uintmax_t maxit = 30;
    //std::cout << "Gamma = " << Gamma << ", cs_sq_0 = " << cs_sq_0 << ", Msq_ratio = " << Msq_ratio << "\n";
    //T result = halley_iterate(functor_2deriv<T>(Gamma, cs_sq_0, Msq_ratio), 0.1, 0.0, 1E+05, get_digits, maxit);
    T result = newton_raphson_iterate(functor_n1deriv<T>(Msq_local, n_ext, cs_to_c_sq_ext, eta, Gamma), 100.0, 0.0, 1E+08, get_digits, maxit);
    //std::cout << "Found res = " << result << "\n";
    return result;
}

#endif //MHD_SOLVER_INCLUDE_EQUATIONS_SOLVERS_H

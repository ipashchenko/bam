#ifndef MHD_SOLVER_INCLUDE_EQUATIONS_SOLVERS_H
#define MHD_SOLVER_INCLUDE_EQUATIONS_SOLVERS_H

#include <iostream>
#include <tuple>
#include <boost/math/tools/roots.hpp>
#include <limits>
// #include "utils.h"
#include <cmath>


// This aimed to solve n_local given M_local, n_ext, c_ext and eta
template <class T>
struct functor_n1deriv
{
        // Functor returning both 1st derivative.
        functor_n1deriv(T const& param_c, T const& param_a, double nu, T const& k_r):
                param_c_(param_c), param_a_(param_a), nu_(nu), k_r_(k_r)
        {}
        std::tuple<T, T> operator()(T const& t)
        {
            // Return both f(x) and f'(x).
            T ft = param_c_ * (asinh(t) + t * sqrt(pow(t, 2) + 1)) - param_a_ * pow(nu_, -1 / k_r_);
            T dt = 2 * param_c_ * sqrt(pow(t, 2) + 1);
            return std::make_tuple(ft, dt);  // 'return' fx, dx.
        }
    private:
        T param_c_;
        T param_a_;
        T nu_;
        T k_r_;
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
T find_n_local_1deriv(T param_c, T param_a, double nu, T k_r)
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
    T result = newton_raphson_iterate(functor_n1deriv<T>(param_c, param_a, nu, k_r), 0.0, -25.0, 25.0, get_digits, maxit);
    std::cout << "Found res = " << result << "\n";
    return result;
}

#endif //MHD_SOLVER_INCLUDE_EQUATIONS_SOLVERS_H

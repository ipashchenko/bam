#include "Data.h"
#include "DNestModel.h"
#include "DNest4.h"
// #include "equations_solvers.h" //optional
#include "Component.h"


using namespace DNest4;

const bool use_logjitter = true;
// const bool image_shifts_from_core = true;


int main(int argc, char** argv)
{
	
    Data::get_instance().load("c1_4.60");
//    Data::get_instance().load("c2_5.00");
    Data::get_instance().load("x1_8.10");
//    Data::get_instance().load("x2_8.43");
    Data::get_instance().load("u1_15.35");
    Data::get_instance().load("k1_23.79");
    // Data::get_instance().load("q1_43.21");
    // Run DNest4
    CommandLineOptions options(argc, argv);
    // Load sampler options from file
    Options sampler_options(options.get_options_file().c_str());
    sampler_options.print(std::cout);
    Sampler<DNestModel> sampler = setup<DNestModel>(options);
    sampler.run();
    return 0;
}

// int main(int argc, char** argv){
//     CoreComponent core(Component::Gaussian);
//     double a = 1, PA = 1, logsize_1 = 1, k_r = 1, k_theta = 1, lognu_max = 1, logS_max = 1, alpha_thick = 1, alpha_thin = 1, p = 1, c = 1;
//     double nu = 4.6;
//     core.set_params(a, PA, logsize_1, k_r, k_theta, lognu_max, logS_max, alpha_thick, alpha_thin, p, c);
//     std::pair<double, double> result = core.get_pos(nu);
//     std::cout << "Ra: " << result.first << ", Dec: " << result.second << std::endl;
//     double result_div = find_n_local_1deriv(a, c, nu, k_r);
//     std::cout << result_div << std::endl;
//     return 0;
// }
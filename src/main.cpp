#include "Data.h"
#include "DNestModel.h"
#include "DNest4.h"
#include "stats.hpp"
#include "gcem.hpp"


using namespace DNest4;


//void test_student() {
//    // parameters
//    double dof = 14.0;
//
//// standard input
//    double dens_val = stats::dt(0.5,dof);
//    double log_dens_val = stats::dt(0.5,dof,true);
//    std::cout << log_dens_val << std::endl;
//}


// On server
//mkdir Release
//cd Release
//cmake -DCMAKE_BUILD_TYPE=Release ..
//make
int main(int argc, char** argv)
{
    Data::get_instance().load("/home/ilya/github/bam/data/good_difmap.txt");
    // set the sampler and run it!
    Sampler<DNestModel> sampler = setup<DNestModel>(argc, argv);
    sampler.run();

    return 0;
}

//int main() {
//    test_student();
//}
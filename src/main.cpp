#include "Data.h"
#include "DNestModel.h"
#include "DNest4.h"

using namespace DNest4;


// On server
//mkdir Release
//cd Release
//cmake -DCMAKE_BUILD_TYPE=Release ..
//make
int main(int argc, char** argv)
{

//    // Run DNest4
//    CommandLineOptions options(argc, argv);
//
//    // Load sampler options from file
//    Options sampler_options(options.get_options_file().c_str());
//    sampler_options.print(std::cout);
////
//    Data::get_instance().load(options.get_data_file());
//    Sampler<DNestModel> sampler = setup<DNestModel>(options);
//    sampler.run();
//
//    //Data::get_instance().load("/home/ilya/github/bam/data/test_60s.txt");
//    Data::get_instance().load("/home/ilya/data/BK150/CBAM/q1/ta_3cg_bam.txt");
//    // set the sampler and run it!
////    Sampler<DNestModel> sampler = setup<DNestModel>(argc, argv);
//    sampler.run();




    // From cbam
    Data::get_instance().load("/home/ilya/data/rfc/CBAM/ta_3cg_bam.txt");
    // Run DNest4
    CommandLineOptions options(argc, argv);
    // Load sampler options from file
    Options sampler_options(options.get_options_file().c_str());
    sampler_options.print(std::cout);
    Sampler<DNestModel> sampler = setup<DNestModel>(options);
    sampler.run();






    return 0;
}
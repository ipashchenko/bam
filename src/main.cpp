#include "Data.h"
#include "DNestModel.h"
#include "DNest4.h"

using namespace DNest4;


int main(int argc, char** argv)
{
	
//	Data::get_instance().load("c1_4.60");
//	Data::get_instance().load("c2_5.00");
//	Data::get_instance().load("x1_8.10");
//	Data::get_instance().load("x2_8.43");
	Data::get_instance().load("u1_15.35");
	Data::get_instance().load("k1_23.79");
	Data::get_instance().load("q1_43.21");
    // Run DNest4
    CommandLineOptions options(argc, argv);
    // Load sampler options from file
    Options sampler_options(options.get_options_file().c_str());
    sampler_options.print(std::cout);
    Sampler<DNestModel> sampler = setup<DNestModel>(options);
    sampler.run();
    return 0;
}
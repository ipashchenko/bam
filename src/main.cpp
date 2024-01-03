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
    Data::get_instance().load("/home/ilya/Downloads/mojave/0851+202/3comp.txt");
    // set the sampler and run it!
    Sampler<DNestModel> sampler = setup<DNestModel>(argc, argv);
    sampler.run();
	
	// Use custom OPTIONS and data file
//	// Run DNest4
//	CommandLineOptions options(argc, argv);
//
//	// Load sampler options from file
//	Options sampler_options(options.get_options_file().c_str());
//	sampler_options.print(std::cout);
//
//	Data::get_instance().load(options.get_data_file());
//	Sampler<DNestModel> sampler = setup<DNestModel>(options);
//	sampler.run();
	
	
    return 0;
}
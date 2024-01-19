#include "Data.h"
#include "DNestModel.h"
#include "DNest4.h"

using namespace DNest4;

const ComponentType component_type = circular;
const bool hyperpriors = false;
const bool use_jitters = false;
const bool use_offsets = false;


int main(int argc, char** argv)
{
//	// OPTIONS from run directory, data file hard coded in main.
//	Data::get_instance().load("/home/ilya/Downloads/mojave/0136+176/2011_07_24/0136+176.u.2011_07_24_60sec_antennas.csv");
//	// set the sampler and run it!
//	Sampler<DNestModel> sampler = setup<DNestModel>(argc, argv);
//	sampler.run();
	
    // Use custom OPTIONS and data file
	// Run DNest4
	CommandLineOptions options(argc, argv);

	// Load sampler options from file
	Options sampler_options(options.get_options_file().c_str());
	sampler_options.print(std::cout);
	
	const std::string& data_file = options.get_data_file();
	Data::get_instance().load(data_file.c_str());
	Sampler<DNestModel> sampler = setup<DNestModel>(options);
	sampler.run();
	
	return 0;
}

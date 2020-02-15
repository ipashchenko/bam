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

    // Run DNest4
    CommandLineOptions options(argc, argv);
    Data::get_instance().load(options.get_data_file());
    Sampler<DNestModel> sampler = setup<DNestModel>(options);
    sampler.run();

    //Data::get_instance().load("/home/ilya/github/bam/data/test_60s.txt");
    ////Data::get_instance().load("/home/ilya/github/bam/data/1502+106.u.2003_03_29.120s.txt");
    //// set the sampler and run it!
    //Sampler<DNestModel> sampler = setup<DNestModel>(argc, argv);
    //sampler.run();

    return 0;
}
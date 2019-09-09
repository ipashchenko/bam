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
    Data::get_instance().load("/home/ilya/github/bsc/test.txt");
    // set the sampler and run it!
    Sampler<DNestModel> sampler = setup<DNestModel>(argc, argv);
    sampler.run();

    return 0;
}
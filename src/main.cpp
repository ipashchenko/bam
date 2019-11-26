#include "Data.h"
#include "DNestModel.h"
#include "DNest4.h"

using namespace DNest4;


// On server
//mkdir Release
//cd Release
//cmake -DCMAKE_BUILD_TYPE=Release ..
//make
// TODO: Check txt-files - some rows with NaNs should be deleted
int main(int argc, char** argv)
{
    Data::get_instance().load("/home/ilya/github/bam/data/1502+106.u.1997_08_18.120s.txt");
    // set the sampler and run it!
    Sampler<DNestModel> sampler = setup<DNestModel>(argc, argv);
    sampler.run();

    return 0;
}
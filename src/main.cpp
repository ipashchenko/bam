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
    Data::get_instance().load("/home/ilya/data/VLBI_Gaia/J1443+0809_C_2013_04_22_pet_vis_noantennas.txt");
    // set the sampler and run it!
    Sampler<DNestModel> sampler = setup<DNestModel>(argc, argv);
    sampler.run();

    return 0;
}

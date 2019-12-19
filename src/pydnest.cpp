#include "Data.h"
#include "DNestModel.h"
#include "DNest4.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using namespace DNest4;
namespace ph = std::placeholders;
namespace py = pybind11;

void run(const char* datafile)
{
    Data::get_instance().load(datafile);
    Options sampler_options(1,
                            50000,
                            50000,
                            100,
                            0,
                            200,
                            100,
                            0);
    Sampler<DNestModel> sampler(5,
                                2.7182818284590451,
                                sampler_options,
                                true,
                                false);

    // Seed RNGs
    uint seed = 32;
    sampler.initialise(seed);

    sampler.run();
}

// Compile with:
// c++ -O3 -Wall -shared -std=c++11 -fPIC `python3 -m pybind11 --includes` -o pydnest`python3-config --extension-suffix` src/pydnest.cpp src/Data.cpp src/Component.cpp src/SkyModel.cpp src/DNestModel.cpp  -I/home/ilya/github/bam/src -I/home/ilya/github/DNest4/code
PYBIND11_MODULE(pydnest, m) {
using namespace pybind11::literals; // for _a literal to define arguments
m.doc() = "DNest for gaussian models"; // optional module docstring

m.def("run", &run, "Sample model using DNest",
"datafile"_a);
}
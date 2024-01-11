#include "Data.h"
#include "DNestModel.h"
#include "DNest4.h"
#include "Gaussian2D.h"

using namespace DNest4;


std::ostream& write_vector(std::ostream &os, std::vector<double> &v) {
	for (auto & i : v) {
		os << i << "\n";
	}
	return os;
}


std::ostream &
write_2dvector(std::ostream &os, std::vector<std::vector<double>> &v) {
	for (auto & i : v) {
		for (double j : i) {
			os << j <<" ";
		}
		os<<"\n";
	}
	return os;
}

// On server
//mkdir Release
//cd Release
//cmake -DCMAKE_BUILD_TYPE=Release ..
//make
int main(int argc, char** argv)
{
    Data::get_instance().load("/home/ilya/Downloads/mojave/0212+735/0212+735.u.2019_08_15_60sec_antennas.txt");
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

int main1(){
	DNest4::RNG rng;
	rng.set_seed(100);
	// MOJAVE
//	Gaussian2D gaussian_2_d({-1.74, -0.71}, -0.41, 1.62, 1.03);
	// Twice narrower
	Gaussian2D gaussian_2_d(-1.74, -0.71, -0.41*2, 1.62, 1.03);
	std::vector<std::vector<double>> samples;
	std::vector<double> sample, cdf, inverse_cdf, u{0.5, 0.5}, x0{3., 0.};
	cdf = gaussian_2_d.cdf(x0);
	std::cout << "cdf[0] = " << cdf[0] << ", cdf[1] = " << cdf[1] << "\n";
	inverse_cdf = gaussian_2_d.cdf_inverse(u);
	std::cout << "inverse_cdf[0] = " << inverse_cdf[0] << ", inverse_cdf[1] = " << inverse_cdf[1] << "\n";

	std::cout << "Before perturb: (" << x0[0] << ", " << x0[1] << ")\n";
	gaussian_2_d.perturb(x0, reinterpret_cast<RNG &>(rng));
	std::cout << "After perturb: (" << x0[0] << ", " << x0[1] << ")\n";

	std::string file("samples.txt");
	std::remove(file.c_str());
	

	for(size_t i = 0; i < 100000; i++){
		sample = gaussian_2_d.generate(rng);
		samples.push_back(sample);
	}
	
	std::fstream fs;
	std::stringstream ss;
	fs.open(file, std::ios::out | std::ios::app);
	if (fs.is_open()) {
		write_2dvector(fs, samples);
		fs.close();
	}
	
	return 0;
}

template< class T, class... Args >
std::shared_ptr<T> make_prior( Args&&... args ) { return std::make_shared<T>(args...); }

int main0(){

//	std::shared_ptr<Gaussian2D> FluxSizePrior;
//	FluxSizePrior = make_prior<Gaussian2D>(0., 0., 0.5, 1., 2.);
	
	std::shared_ptr<Gaussian2D> FluxSizePrior;
	if(!FluxSizePrior) {
		FluxSizePrior = std::make_shared<Gaussian2D>(0., 0., 0.5, 1., 2.);
	}
	auto res = FluxSizePrior->cdf({1., 0.5});
	std::cout << res[0] << ", " << res[1] << "\n";
	
}
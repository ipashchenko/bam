#include "Data.h"
#include <fstream>
#include <vector>
#include <iostream>
#include <set>
#include <iomanip>
#include "rapidcsv.h"


using namespace std;

// The static instance
Data Data::instance;

Data::Data() = default;


void Data::load(const char* filename)
{
	
	rapidcsv::Document doc(filename);
	ant_i = doc.GetColumn<int>("t1");
	ant_j = doc.GetColumn<int>("t2");
	std::vector<double> _u = doc.GetColumn<double>("u");
	std::vector<double> _v = doc.GetColumn<double>("v");
	std::vector<double> _vis_real = doc.GetColumn<double>("vis_re");
	std::vector<double> _vis_imag = doc.GetColumn<double>("vis_im");
	std::vector<double> _sigma = doc.GetColumn<double>("error");

	std::cout << "First raw : " << std::setprecision(15) << ant_i[0] << ", " << ant_j[0] << ", " << _u[0] << ", " << _v[0] << ", " << _vis_real[0] << ", " << _vis_imag[0] << ", " << _sigma[0] << "\n";
	
	std::vector<int> _antennas;
	
	// Vector of the unique antenna numbers
	antennas.insert(antennas.end(), ant_i.begin(), ant_i.end());
	antennas.insert(antennas.end(), ant_j.begin(), ant_j.end());
	std::set<int> s(antennas.begin(), antennas.end());
	antennas.assign(s.begin(), s.end());
	
	// Generate the map between ant_i and its position in antennas vector
	for (int i=0; i<antennas.size();i++) {
		antennas_map[antennas[i]] = i;
	}
	
	std::cout << "Antenna map :" << std::endl;
	for (auto x : antennas_map) {
		std::cout << x.first << " -- " << x.second << std::endl;
	}
	
	// Generate the map between antenna position in antennas vector and ant_i
	for (int i=0; i<antennas.size();i++) {
		antennas_map_inv[i] = antennas[i];
	}
	
	std::cout << "Antenna map inverse :" << std::endl;
	for (auto x : antennas_map_inv) {
		std::cout << x.first << " -- " << x.second << std::endl;
	}
	
    // Copy the data to the valarrays
    u = valarray<double>(&_u[0], _u.size());
    v = valarray<double>(&_v[0], _v.size());
    vis_real = valarray<double>(&_vis_real[0], _vis_real.size());
    vis_imag = valarray<double>(&_vis_imag[0], _vis_imag.size());
    sigma = valarray<double>(&_sigma[0], _sigma.size());
}

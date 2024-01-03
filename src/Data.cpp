#include "Data.h"
#include <fstream>
#include <vector>
#include <iostream>
#include <set>

using namespace std;

// The static instance
Data Data::instance;

Data::Data() = default;

// TODO: Use files with header
void Data::load(const char* filename)
{
	std::vector<int> _antennas;
    // Vectors to hold the data
    std::vector<double> _u;
    std::vector<double> _v;
    std::vector<double> _vis_real;
    std::vector<double> _vis_imag;
    std::vector<double> _sigma;

    // Open the file
    fstream fin(filename, ios::in);

    // Temporary variables
    // u, v, re, im, var
    double temp1, temp2, temp3, temp4, temp5;
	int temp01, temp02;

    // Read until end of file
    while(fin >> temp01 && fin >> temp02 && fin >> temp1 && fin >> temp2 && fin >> temp3 && fin >> temp4 && fin >> temp5)
    {
		ant_i.push_back(temp01);
		ant_j.push_back(temp02);
        _u.push_back(temp1);
        _v.push_back(temp2);
        _vis_real.push_back(temp3);
        _vis_imag.push_back(temp4);
        _sigma.push_back(temp5);
    }

    // Close the file
    fin.close();
	
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

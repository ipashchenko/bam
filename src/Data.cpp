#include "Data.h"
#include "utils.h"
#include <fstream>
#include <vector>
#include <iostream>
#include <set>

using namespace std;

// The static instance
Data Data::instance;

Data::Data() = default;

void Data::load(const std::string& filename)
{
	
	std::vector<std::string> splitted = split(filename, "_");
	std::string band, frequency_ghz;
	band = splitted[0];
	frequency_ghz = splitted[1];
	bands.push_back(band);
	band_freq_map[band] = std::stod(frequency_ghz);
	
	cout << "Loading file " << filename << ", band : " << band << ", Freq[GHz] = " << frequency_ghz << endl;
	
    // Vectors to hold the data
    std::vector<double> _u;
    std::vector<double> _v;
	std::vector<std::complex<double>> _vis;
    std::vector<double> _vis_real;
    std::vector<double> _vis_imag;
    std::vector<double> _sigma;

    // Open the file
    fstream fin(filename, ios::in);

    // Temporary variables
    // u, v, re, im, sigma
    double temp1, temp2, temp3, temp4, temp5;

    // Read until end of file
    while(fin >> temp1 && fin >> temp2 && fin >> temp3 && fin >> temp4 && fin >> temp5)
    {
        _u.push_back(temp1);
        _v.push_back(temp2);
        _vis_real.push_back(temp3);
        _vis_imag.push_back(temp4);
		_vis.push_back({temp3, temp4});
        _sigma.push_back(temp5);
    }

    // Close the file
    fin.close();

	cout << "Loaded " << _u.size() << " visibilities\n";
	
    // Copy the data to the containers
	Eigen::Map<ArrayXd> u_(_u.data(), _u.size());
	Eigen::Map<ArrayXd> v_(_v.data(), _v.size());
    u[band] = u_;
    v[band] = v_;
	Eigen::Map<ArrayXcd> vis_(_vis.data(), _vis.size());
    vis[band] = vis_;
	Eigen::Map<ArrayXd> sigma_(_sigma.data(), _sigma.size());
    sigma[band] = sigma_;
}

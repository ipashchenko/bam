#ifndef BAM_DATA_H
#define BAM_DATA_H

#include <vector>
#include <set>
#include <unordered_map>
#include <string>
#include <eigen3/Eigen/Core>

using Eigen::ArrayXcd;
using Eigen::ArrayXd;

class Data
{
    private:
        // Static "global" instance
        static Data instance;
		std::vector<std::string> bands;
		std::unordered_map<std::string, double> band_freq_map;
        // uv-coordinates
        std::unordered_map<std::string, ArrayXd> u;
        std::unordered_map<std::string, ArrayXd> v;
        // visibilities
        std::unordered_map<std::string, ArrayXcd> vis;
        // error of real/imag part
        std::unordered_map<std::string, ArrayXd> sigma;

    public:
        // Constructor
        Data();
        // Getter for the global instance
        static Data& get_instance()
        { return instance; }

        // Load data from a file
        void load(const std::string& filename);

        // Access to the data points
        ArrayXd& get_u(std::string band)
        { return u[band]; }
        ArrayXd& get_v(std::string band)
        { return v[band]; }
        ArrayXcd& get_vis_real(std::string band)
        { return vis[band]; }
        ArrayXd& get_sigma(std::string band)
        { return sigma[band]; }
		std::unordered_map<std::string, double> get_band_freq_map()
		{ return band_freq_map; }
		std::vector<std::string> get_bands()
		{ return bands; }
};


#endif //BAM_DATA_H

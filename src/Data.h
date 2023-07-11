#ifndef BAM_DATA_H
#define BAM_DATA_H

#include <valarray>
#include <vector>
#include <set>
#include <unordered_map>
#include <string>

class Data
{
    private:
        // Static "global" instance
        static Data instance;
		std::vector<std::string> bands;
		std::unordered_map<std::string, double> band_freq_map;
        // uv-coordinates
        std::unordered_map<std::string, std::valarray<double>> u;
        std::unordered_map<std::string, std::valarray<double>> v;
        // visibilities
        std::unordered_map<std::string, std::valarray<double>> vis_real;
        std::unordered_map<std::string, std::valarray<double>> vis_imag;
        // error of real/imag part
        std::unordered_map<std::string, std::valarray<double>> sigma;

    public:
        // Constructor
        Data();
        // Getter for the global instance
        static Data& get_instance()
        { return instance; }

        // Load data from a file
        void load(const std::string& filename);

        // Access to the data points
        std::valarray<double>& get_u(std::string band)
        { return u[band]; }
        std::valarray<double>& get_v(std::string band)
        { return v[band]; }
        std::valarray<double>& get_vis_real(std::string band)
        { return vis_real[band]; }
        std::valarray<double>& get_vis_imag(std::string band)
        { return vis_imag[band]; }
        std::valarray<double>& get_sigma(std::string band)
        { return sigma[band]; }
		std::unordered_map<std::string, double> get_band_freq_map()
		{ return band_freq_map; }
		std::vector<std::string> get_bands()
		{ return bands; }
};


#endif //BAM_DATA_H

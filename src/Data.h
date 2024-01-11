#ifndef BAM_DATA_H
#define BAM_DATA_H

#include <valarray>
#include <vector>
#include <set>
#include <unordered_map>


class Data
{
    private:
        // Static "global" instance
        static Data instance;
		// Unique antenna numbers
		std::vector<int> antennas;
		// Map between ant_i and its position in antennas vector
		std::unordered_map<int, int> antennas_map;
		// Map between antenna position in antennas vector and ant_i
		std::unordered_map<int, int> antennas_map_inv;
		// Antenna numbers
		std::vector<int> ant_i;
		std::vector<int> ant_j;
        // uv-coordinates
        std::valarray<double> u;
        std::valarray<double> v;
        // visibilities
        std::valarray<double> vis_real;
        std::valarray<double> vis_imag;
        // error of real/imag part
        std::valarray<double> sigma;

    public:
        // Constructor
        Data();
        // Getter for the global instance
        static Data& get_instance()
        { return instance; }
		
		// Number of antennas
		int n_antennas() const
		{ return antennas.size(); }
		
		int n_vis() const
		{ return u.size(); }
		
        // Load data from a file
        void load(const char* filename);

        // Access to the data points
		const std::vector<int>& get_ant_i() const
		{ return ant_i; }
		const std::vector<int>& get_ant_j() const
		{ return ant_j; }
        const std::valarray<double>& get_u() const
        { return u; }
        const std::valarray<double>& get_v() const
        { return v; }
        const std::valarray<double>& get_vis_real() const
        { return vis_real; }
        const std::valarray<double>& get_vis_imag() const
        { return vis_imag; }
        const std::valarray<double>& get_sigma() const
        { return sigma; }
		std::unordered_map<int, int>& get_antennas_map()
		{ return antennas_map; }
		std::unordered_map<int, int>& get_antennas_map_inv()
		{ return antennas_map_inv; }
};


#endif //BAM_DATA_H

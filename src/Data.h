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

        // Load data from a file
        void load(const std::string& filename);

        // Access to the data points
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
};


#endif //BAM_DATA_H

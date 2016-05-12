#include <string>
#include <fstream>

using namespace std;

class FluxDensityInfo {
    private:
        int num_flux_;
        double *read_tmp;        
        void dealloc();

    public:
        //double *read_tmp;
        double flux_min, flux_max;

        double *S, *S2p5dNdS, *dNdS, *logS;
        double *dNdlogS;

        FluxDensityInfo (string &FileName, double flux_min_mJy, double flux_max_mJy);
        ~FluxDensityInfo();

        /* returns the number of flux */
        int Num_Flux() { return num_flux_; }

        void Set(int nflux, int ncol);
        void Read(ifstream &infile);
};

class rngHandle;
class FluxDensityInfo;
template<typename T> class Healpix_Map;

template<typename T> void create_poisson_map (rngHandle &rng, FluxDensityInfo &flux, Healpix_Map<T> &map);

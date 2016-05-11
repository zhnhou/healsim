#include <fstream>

class FluxDensityInfo {
    private:
        int num_flux_;
        double *read_tmp;        
        void dealloc();

    public:
        //double *read_tmp;
        int    nside, npix;
        double flux_min, flux_max;

        double *S, *S2p5dNdS, *dNdS, *logS;
        double *dNdlogS;

        FluxDensityInfo (string &FileName);
        ~FluxDensityInfo();

        /* returns the number of flux */
        int Num_Flux() { return num_flux_; }

        void Set(int nflux, int ncol);

        void Read(ifstream &infile);
};

#include <fstream>

using namespace std;

class FluxDensity {
    private:
        int num_flux_;
        double *read_tmp;        
        void dealloc();

    public:
        //double *read_tmp;
        double *S, *S2p5dNdS, *dNdS, *logS;
        double *dNdlogS;

        FluxDensity (string &FileName);
        ~FluxDensity();

        /* returns the number of flux */
        int Num_Flux() { return num_flux_; }

        void Set(int nflux, int ncol);

        void Read(ifstream &infile);
};

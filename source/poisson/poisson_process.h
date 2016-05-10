#include "arr.h"

using namespace std;

class FluxDensity {
    private:
        arr<double> read_tmp_, S_, S2p5dNdS_, dNdS_, logS_, dNdlogS_;
        int num_flux;
        void dealloc();

    public:
        FluxDensity (string &FileName);
        ~FluxDensity();

        void assertArraySizes() const;

        /* returns the number of flux */
        int Num_Flux() { return num_flux; }
        /* returns the S array */
        arr<double> &S() { return S_; }
        /* returns the $S^2.5 dN/dS$ array */
        arr<double> &S2p5dNdS() { return S2p5dNdS_; }
        /* returns the dN/dS array */
        arr<double> &dNdS() { return dNdS_; }
        /* returns the logS array */
        arr<double> &logS() { return logS_; }
        /* returns the dN/dlogS array */
        arr<double> &dNdlogS() { return dNdlogS_; }

        void Set(int nflux);
};

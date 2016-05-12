#include <string>
#include "mkl_vsl.h"

using namespace std;

class rngHandle {
    private:
        string  rng_path_;
        int     seed_;
        bool    mkl_rng_, hpx_rng_;
    public:
        rngHandle (int seed, string &rng_cache_path, bool mkl_rng, bool hpx_rng);
        ~rngHandle();

        VSLStreamStatePtr vsl_poisson_stream;
        VSLStreamStatePtr vsl_uniform_stream;
        
        bool vsl_poisson_init;
        bool vsl_uniform_init;
        bool hpx_gaussian_init;

        int  Set(int isim);
        void save_mkl_rng (string &rng_file);
        void read_mkl_rng (string &rng_file);
};

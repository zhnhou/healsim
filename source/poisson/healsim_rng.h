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

        VSLStreamStatePtr vsl_stream_poisson;
        VSLStreamStatePtr vsl_stream_uniform;
        
        int     vsl_bng;
        int     vsl_method_poisson, vsl_method_uniform;

        bool vsl_poisson_init;
        bool vsl_uniform_init;
        bool hpx_gaussian_init;

        int  Set(int isim);
        void save_mkl_rng (string &rng_file);
        void read_mkl_rng (string &rng_file);
};

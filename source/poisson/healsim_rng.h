#include <string>
#include "mkl_vsl.h"

using namespace std;

class rngHandle {
    private:
        string rng_path_;
    public:
        rngHandle (int seed, int isim, string &rng_cache_path, bool mkl_rng, bool hpx_rng);
        ~rngHandle();

        VSLStreamStatePtr vsl_poisson_stream;
        VSLStreamStatePtr vsl_uniform_stream;
        
        bool vsl_poisson_init;
        bool vsl_uniform_init;
        bool hpx_rng_init;

        void save_mkl_rng (string &rng_file);
        void read_mkl_rng (string &rng_file);
};

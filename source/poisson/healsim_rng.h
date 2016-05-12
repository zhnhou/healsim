#include <string>

using namespace std;

class rngHandle {
    private:
        string rng_path_;
    public:
        rngHandle (int seed, int isim, string &rng_cache_path, bool mkl_rng, bool hpx_rng);
        ~rngHandle();

        VSLStreamStatePtr vsl_poisson_stream;
        VSLStreamStatePtr vsl_uniform_stream;
        
        bool vsl_poisson_init=false;
        bool vsl_uniform_init=false;
        bool hpx_rng_init=false;

        void save_mkl_rng (string &rng_file);
        void read_mkl_rng (string &rng_file);
};

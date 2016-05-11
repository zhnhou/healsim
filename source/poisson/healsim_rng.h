#include <string>

class rngHandle {
    private:
        std::string rng_path_;
    public:
        VSLStreamStatePtr vsl_poisson_stream;
        VSLStreamStatePtr vsl_uniform_stream;
        
        bool vsl_poisson_init=false;
        bool vsl_uniform_init=false;
        bool hpx_rng_init=false;

        rngHandle (int seed, int isim, std::string &rng_cache_path, bool mkl_rng, bool hpx_rng);
        ~rngHandle();

        void save_mkl_rng (std::string &rng_file);
        void retrieve ();
};

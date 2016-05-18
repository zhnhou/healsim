#include "paramfile.h"
#include "announce.h"
#include "healpix_map.h"

#include "mkl_vsl.h"
#include "healsim_rng.h"
#include "poisson_process.h"

using namespace std;

namespace {
template<typename T> void syn_poisson (int isim, rngHandle &rng, FluxDensityInfo &flux, Healpix_Map<T> &poisson_map, int feedback=0) {

    rng.Set(isim);
    create_poisson_map (rng, flux, poisson_map);

}
} // unnamed namespace

int syn_poisson_module (int argc, const char **argv) {
    module_startup ("syn_poisson", argc, argv);
    paramfile params (getParamsFromCmdline(argc, argv));

    int  nsim = params.find<int> ("num_sim", 1);
    int  feedback = params.find<int> ("feedback", 0);
    bool dp = params.find<bool> ("double_precision", false);

    int nside = params.find<int>("nside");
    int seed  = params.find<int>("rand_seed");
    string infile = params.find<string>("input_file");
    string outpath = params.find<string>("output_path");
    string outpre  = params.find<string>("output_prefix");
    
    double flux_min = params.find<double>("flux_min_mJy", 0.00e0);
    double flux_max = params.find<double>("flux_max_mJy");

    string rng_cache_path = outpath + "/rng_cache";
    rngHandle rng (seed, rng_cache_path, true, false);
    if (feedback > 0)
        cout << "MKL RNG initialized" << endl;

    FluxDensityInfo flux (infile, flux_min, flux_max);
    if (feedback > 0)
        cout << "Flux Density initialized" << endl;

    if (dp) {
        Healpix_Map<double> poisson_map(nside, RING, SET_NSIDE);
        for (int isim=0; isim<nsim; isim++)
            syn_poisson<double>(isim, rng, flux, poisson_map, feedback);}
    else {
        Healpix_Map<float> poisson_map(nside, RING, SET_NSIDE);
        for (int isim=0; isim<nsim; isim++)
            syn_poisson<float>(isim, rng, flux, poisson_map, feedback);}

    return 0;
}

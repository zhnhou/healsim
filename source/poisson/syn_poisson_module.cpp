#include "paramfile.h"
#include "announce.h"
#include "healpix_map.h"

#include "mkl_vsl.h"
#include "healsim_rng.h"
#include "poisson_process.h"

using namespace std;

namespace {
template<typename T> void syn_poisson (paramfile &params, int isim, int feedback=0) {

    int nside = params.find<int>("nside");
    int seed  = params.find<int>("rand_seed");
    string infile = params.find<string>("input_file");
    string outpath = params.find<string>("output_path");
    string outpre  = params.find<string>("output_prefix");
    
    double flux_min = params.find<double>("flux_min_mJy", 0.00e0);
    double flux_max = params.find<double>("flux_max_mJy");

    string rng_cache_path = outpath + "/rng_cache";

    rngHandle rng (seed, isim, rng_cache_path, true, false);
    FluxDensityInfo flux (infile);
    flux.flux_min = flux_min;
    flux.flux_max = flux_max;

    Healpix_Map<T> map(nside, RING, SET_NSIDE);

    if (feedback > 0)
        cout << "Flux Density initialized" << endl;
    
    //cout << flux.S()[0] << endl;
}
} // unnamed namespace

int syn_poisson_module (int argc, const char **argv) {
    module_startup ("syn_poisson", argc, argv);
    paramfile params (getParamsFromCmdline(argc, argv));

    int  nsim = params.find<int> ("num_sim", 1);
    int  feedback = params.find<int> ("feedback", 0);
    bool dp = params.find<bool> ("double_precision", false);

    for (int isim=0; isim<nsim; isim++) {
        dp ? syn_poisson<double>(params, isim, feedback) : syn_poisson<float>(params, isim, feedback);
    }
    return 0;
}

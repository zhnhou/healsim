#include "paramfile.h"
#include "announce.h"
#include "poisson_process.h"

using namespace std;

namespace {
template<typename T> void syn_poisson (paramfile &params, int isim) {

    int nside = params.find<int>("nside");
    int seed  = params.find<int>("rand_seed");
    string infile = params.find<string>("input_file");
    string outpath = params.find<string>("output_path");
    string outpre  = params.find<string>("output_prefix");
    
    double flux_min = params.find<double>("flux_min_mJy", 0.00e0);
    double flux_max = params.find<double>("flux_max_mJy");

    FluxDensity flux (infile);
    
    //cout << flux.S()[0] << endl;
}
} // unnamed namespace

int syn_poisson_module (int argc, const char **argv) {
    module_startup ("syn_poisson", argc, argv);
    paramfile params (getParamsFromCmdline(argc, argv));

    int  nsim = params.find<int> ("num_sim", 1);
    bool dp = params.find<bool> ("double_precision", false);

    for (int isim=0; isim<nsim; isim++) {
        dp ? syn_poisson<double>(params, isim) : syn_poisson<float>(params, isim);
    }
    return 0;
}

#include "paramfile.h"
#include "announce.h"

using namespace std;

namespace {

template<typename T> void syn_poisson (paramfile &params) {

    int nside = params.template find<int>("nside");
    int seed  = params.template find<int>("rand_seed");
    string infile = params.template find<string>("input_file");
    string outpath = params.template find<string>("output_path");
    string outpre  = params.template find<string>("output_prefix");
    
    double flux_min = params.template find<double>("flux_min_mJy", 0.00e0);
    double flux_max = params.template find<double>("flux_max_mJy");
}
} // unnamed namespace

int syn_poisson_module (int argc, const char **argv) {
    module_startup ("syn_poisson", argc, argv);
    paramfile params (getParamsFromCmdline(argc, argv));

    bool dp = params.find<bool> ("double_precision",false);
    dp ? syn_poisson<double>(params) : syn_poisson<float>(params);
    return 0;
}

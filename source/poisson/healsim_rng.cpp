#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "mkl_vsl.h"
#include "healsim_rng.h"

rngHandle::rngHandle(int seed, int isim, string &rng_cache_path, bool mkl_rng, bool hpx_rng) {
    rng_path_ = rng_cache_path;

    if (hpx_rng) {


    } else if (mkl_rng) {
        stringstream ss;
        ss << isim;
        string vsl_poisson_rng = rng_path_+"/vsl_poisson_sim_"+ss.str()+".rng";

        ifstream in_rng (vsl_poisson_rng);

        if (in_rng.good()) {
            retrieve();
        } else if ( (! vsl_poisson_init) || (! vsl_uniform_init) ) {
            vslNewStream( &vsl_poisson_stream, VSL_BRNG_MCG31, seed+10*isim);
            vslNewStream( &vsl_uniform_stream, VSL_BRNG_MCG31, -1*(seed-100*isim));
            save_mkl_rng(vsl_poisson_rng);
        }

        vsl_poisson_init = true;
        vsl_uniform_init = true;
    } else {

    }
}

rngHandle::~rngHandle() {

}

void rngHandle::save_mkl_rng(string &rng_file) {
    ofstream outbin(rng_file.c_str(), ios::out | ios::binary);
    outbin.write( (char*) &vsl_poisson_stream, sizeof(vsl_poisson_stream) );
    outbin.write( (char*) &vsl_uniform_stream, sizeof(vsl_uniform_stream) );
    outbin.close();
}

void rngHandle::retrieve() {

}

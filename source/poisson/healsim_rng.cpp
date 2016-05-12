#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>

#include "mkl_vsl.h"
#include "healsim_rng.h"

rngHandle::rngHandle(int seed, string &rng_cache_path, bool mkl_rng, bool hpx_rng) {
    rng_path_ = rng_cache_path;
    seed_ = seed;
    mkl_rng_ = mkl_rng;
    hpx_rng_ = hpx_rng;

    if (mkl_rng && hpx_rng) {
        cout << "Both mkl_rng and hpx_rng are true. STOP!" << endl;
        exit(1);
    } else if (! (mkl_rng && hpx_rng)) {
        cout << "Neither mkl_rng or hpx_rng is true. STOP!" << endl;
        exit(1);
    }

    vsl_poisson_init = false;
    vsl_uniform_init = false;
    hpx_gaussian_init = false;
}

int rngHandle::Set(int isim) {
    if (hpx_rng_) {


    } else if (mkl_rng_) {
        stringstream ss;
        ss << isim;
        string vsl_poisson_rng = rng_path_+"/vsl_poisson_sim_"+ss.str()+".rng";

        ifstream in_rng (vsl_poisson_rng.c_str());

        if (in_rng.good()) {
            read_mkl_rng(vsl_poisson_rng);
        } else if ( (! vsl_poisson_init) || (! vsl_uniform_init) ) {
            vslNewStream( &vsl_poisson_stream, VSL_BRNG_MCG31, seed_+10*isim);
            vslNewStream( &vsl_uniform_stream, VSL_BRNG_MCG31, -1*(seed_-100*isim));
            save_mkl_rng(vsl_poisson_rng);
        }

        vsl_poisson_init = true;
        vsl_uniform_init = true;

        return 0;
    }

}

rngHandle::~rngHandle() {

}

void rngHandle::save_mkl_rng(string &rng_file) {
    ofstream outbin(rng_file.c_str(), ios::out | ios::binary);
    if (outbin.good()) {
        outbin.write( (char*) &vsl_poisson_stream, sizeof(vsl_poisson_stream) );
        outbin.write( (char*) &vsl_uniform_stream, sizeof(vsl_uniform_stream) );
    } else {
        cout << "Error in writting to file - " << endl;
        cout << rng_file << endl;
        exit(1);
    }
    outbin.close();
}

void rngHandle::read_mkl_rng(string &rng_file) {
    ifstream inbin(rng_file.c_str(), ios::in | ios::binary);

    if (inbin.good()) {
        inbin.read( (char*) &vsl_poisson_stream, sizeof(vsl_poisson_stream) );
        inbin.read( (char*) &vsl_uniform_stream, sizeof(vsl_uniform_stream) );
    } else {
        cout << "Error in reading from file - " << endl;
        cout << rng_file << endl;
        exit(1);
    }
    inbin.close();
}

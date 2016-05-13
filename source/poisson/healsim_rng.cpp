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
    } else if (! (mkl_rng || hpx_rng)) {
        cout << "Neither mkl_rng or hpx_rng is true. STOP!" << endl;
        exit(1);
    }

    vsl_bng = VSL_BRNG_MCG31;
    vsl_method_poisson = VSL_RNG_METHOD_POISSONV_POISNORM;
    vsl_method_uniform = VSL_RNG_METHOD_UNIFORM_STD;

    vsl_poisson_init = false;
    vsl_uniform_init = false;
    hpx_gaussian_init = false;
}

int rngHandle::Set(int isim) {
    if (hpx_rng_) {


    } else if (mkl_rng_) {
        VSLBRngProperties properties;
        vslGetBrngProperties(VSL_BRNG_MCG31, &properties);
        vsl_stream_size_ = properties.StreamStateSize;

        stringstream ss;
        ss << isim;
        string vsl_poisson_rng = rng_path_+"/vsl_poisson_sim_"+ss.str()+".rng";

        ifstream in_rng (vsl_poisson_rng.c_str());

        if (in_rng.good()) {
            read_mkl_rng(vsl_poisson_rng);
        } else if ( (! vsl_poisson_init) || (! vsl_uniform_init) ) {
            vslNewStream( &vsl_stream_poisson, VSL_BRNG_MCG31, seed_+10*isim+100);
            vslNewStream( &vsl_stream_uniform, VSL_BRNG_MCG31, -1*(seed_-100*isim));
            save_mkl_rng(vsl_poisson_rng);
        }

        vsl_poisson_init = true;
        vsl_uniform_init = true;

        return 0;
    }

}

rngHandle::~rngHandle() {

}

/*
For saving and retrieving the vsl stream, here is a very good direction
https://software.intel.com/en-us/forums/intel-math-kernel-library/topic/309910
*/

void rngHandle::save_mkl_rng(string &rng_file) {
    ofstream outbin(rng_file.c_str(), ios::out | ios::binary);
    if (outbin.good()) {
        outbin.write( (char*) vsl_stream_poisson, vsl_stream_size_ );
        outbin.write( (char*) vsl_stream_uniform, vsl_stream_size_ );
    } else {
        cout << "Error in writting to file - " << endl;
        cout << rng_file << endl;
        exit(1);
    }
    outbin.close();
}

void rngHandle::read_mkl_rng(string &rng_file) {
    ifstream inbin(rng_file.c_str(), ios::in | ios::binary);

    VSLStreamStatePtr vsl_stream_poisson_cache;
    vslNewStream( &vsl_stream_poisson_cache, VSL_BRNG_MCG31, 0);

    VSLStreamStatePtr vsl_stream_uniform_cache;
    vslNewStream( &vsl_stream_uniform_cache, VSL_BRNG_MCG31, 0);

    if (inbin.good()) {
        /* Pay special attention here:
           it is (char*) vsl_stream_poisson_cache, instead of (char*) &vsl_stream_poisson_cache
           a relevant explanation can be found at
           http://www.cplusplus.com/forum/general/76997/

           same treatment in outbin.write()
        */
        inbin.read( (char*) vsl_stream_poisson_cache, vsl_stream_size_ );
        inbin.read( (char*) vsl_stream_uniform_cache, vsl_stream_size_ );

        vslNewStream( &vsl_stream_poisson, VSL_BRNG_MCG31, 0);
        vslNewStream( &vsl_stream_uniform, VSL_BRNG_MCG31, 0);

        //status = vslCopyStreamState( deststream, srcstream );
    
        vslCopyStreamState( vsl_stream_poisson, vsl_stream_poisson_cache );
        vslCopyStreamState( vsl_stream_uniform, vsl_stream_uniform_cache );
    } else {
        cout << "Error in reading from file - " << endl;
        cout << rng_file << endl;
        exit(1);
    }
    inbin.close();
}

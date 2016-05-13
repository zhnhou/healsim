#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <time.h>
#include <sys/time.h>

#include "lsconstants.h"
#include "healpix_map.h"
#include "healsim_rng.h"
#include "poisson_process.h"


FluxDensityInfo::FluxDensityInfo(string &FileName, double flux_min_mJy=0.00e0, double flux_max_mJy=50.00e0) {

    string      line;
    ifstream    infile(FileName.c_str());

    flux_min = flux_min_mJy;
    flux_max = flux_max_mJy;

    if (not infile.is_open())
        cout << "File not open\n";

    int ct=0;
    int num_col=0, num_row=0;
    double n;
    while (! infile.eof()) {
        ct = 0;
        getline(infile, line);

        stringstream icol(line);
        while (icol >> n) ct++;

        if (num_row == 0) num_col = ct;
        if (num_col == ct) num_row++;
    }
    infile.close();

    Set(num_row, num_col);

    infile.open(FileName.c_str());
    Read(infile);
    infile.close();
    
}

FluxDensityInfo::~FluxDensityInfo(){}


void FluxDensityInfo::Set(int nflux, int ncol) {
    num_flux_ = nflux;
    S = new double [num_flux_];
    S2p5dNdS = new double [num_flux_];
    dNdS = new double [num_flux_];
    logS = new double [num_flux_];
    dNdlogS = new double [nflux];
    read_tmp = new double [ncol];
    dlogS = S[1] - S[0];
}

void FluxDensityInfo::Read(ifstream& infile) {

    int    ct=0, irow=0;
    double tmp;
    string line;

    while (!infile.eof()) {
        ct = 0;
        getline(infile, line);

        // Here we check if there is empty line, actually this will avoid
        // the loop continue after reading the end of the file. Otherwise,
        // irow will exceed the array bundary causing the following error
        // message:
        //*** glibc detected *** ./syn_poisson: double free or corruption (out)

        if (line.empty()) continue;
        
        stringstream icol(line);
        while (icol >> tmp) {
            read_tmp[ct] = tmp;
            ct++;

            //if we want to do some formatted output,here is an example
            //Note that the format set up should keep within the loop
            //cout.precision(9);
            //cout.width(18);
            //cout << right << scientific << tmp;
            //cin.get();
            
        }
        S[irow] = pow(10.00e0, read_tmp[0]);
        S2p5dNdS[irow] = pow(10.00e0, read_tmp[1]);
        dNdS[irow] = S2p5dNdS[irow] / pow(S[irow], 2.50e0);
        logS[irow] = log(S[irow]);
        dNdlogS[irow] = dNdS[irow] * S[irow];

        irow++;
    }
}

/*
    double wall0 = get_wall_time();
    double cpu0 = get_cpu_time();

    double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}
    double get_cpu_time(){
    return (double)clock() / CLOCKS_PER_SEC;
}

    double wall0 = get_wall_time();
    double cpu0 = get_cpu_time();

    double wall1 = get_wall_time();
    double cpu1 = get_cpu_time();

    cout << "Wall Time = " << wall1 - wall0 << endl;
    cout << "CPU Time  = " << cpu1  - cpu0  << endl;

*/

template<typename T> void create_poisson_map(rngHandle &rng, FluxDensityInfo &flux, Healpix_Map<T> &map) {

    double frac, frac0;
    double S_max = max(flux.flux_max, flux.S[flux.Num_Flux()-1]);

#pragma omp parallel
{
    int ipix;
#pragma omp for schedule (dynamic,5000)
    for (ipix=0; ipix<map.Npix(); ++ipix) map[ipix] = 0.00e0;
}

    int ct=0, ct0=0;
    for (int i=0; i<flux.Num_Flux(); ++i) {
        if (flux.S[i] <= flux.flux_max) ++ct;
        if (flux.S[i] <  flux.flux_min) ++ct0;
    }

    double *dN  = new double [ct+1]();
    double *dN0 = new double [ct0+1]();
    int *dN_Poisson = new int [ct+1]();

    frac = (log(S_max) - log(flux.S[ct-1])) / (log(flux.S[ct]) - log(flux.S[ct-1]));
    if (ct0 != 0) 
        frac0 = (log(flux.flux_min) - log(flux.S[ct0-1])) / (log(flux.S[ct0]) - log(flux.S[ct0-1]));

    for (int i=0; i<ct; ++i) dN[i] = flux.dNdlogS[i] * flux.dlogS * fourpi;
    dN[ct] = flux.dNdlogS[ct] * flux.dlogS * fourpi * frac;

    if (ct0 == 1) {
        dN0[0] = flux.dNdlogS[0] * flux.dlogS * fourpi * frac0;
        dN[0] -= dN0[0]; }
    else if (ct0 > 1) {
        for (int i=0; i<ct0; ++i) dN0[i] = flux.dNdlogS[i] * flux.dlogS * fourpi;
        dN0[ct0] = flux.dNdlogS[ct0] * flux.dlogS * fourpi * frac0;

        for (int i=0; i<=ct; ++i) dN[i] -= dN0[i];
    }

    //int errcode = viRngPoissonV(rng.vsl_method_poisson, rng.vsl_stream_poisson, ct-ct0+1, &dN_Poisson[ct0], &dN[ct0]);
    int errcode = viRngPoissonV(rng.vsl_method_poisson, rng.vsl_stream_poisson, ct+1, dN_Poisson, dN);

}

template void create_poisson_map(rngHandle &rng, FluxDensityInfo &flux, Healpix_Map<double> &map);
template void create_poisson_map(rngHandle &rng, FluxDensityInfo &flux, Healpix_Map<float> &map);

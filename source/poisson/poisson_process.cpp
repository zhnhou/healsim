#include <iostream>
#include <algorithm>
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

using namespace std;


FluxDensityInfo::FluxDensityInfo(string &FileName, double flux_min_mJy=0.00e0, double flux_max_mJy=50.00e0) {

    string      line;
    ifstream    infile(FileName.c_str());

    flux_min = flux_min_mJy * 1.00e-3; // convert to 
    flux_max = flux_max_mJy * 1.00e-3;

    if (not infile.is_open()) {
        cout << "File not open\n";
        exit(1);
    }

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

    dlogS = logS[1] - logS[0];
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

PoissonN::PoissonN(rngHandle &rng, FluxDensityInfo &flux, int npix) {
    Set(flux, npix);
    GenPN(rng);
}

PoissonN::~PoissonN() {}

void PoissonN::Set(FluxDensityInfo &flux, int npix) {

    if (! init) {
        ct_ = 0;
        ct0_ = 0;

        for (int i=0; i<flux.Num_Flux(); ++i) {
            if (flux.S[i] <= flux.flux_max) ++ct_;
            if (flux.S[i] <  flux.flux_min) ++ct0_;
        }

        dN  = new double [ct_+1]();
        dN0 = new double [ct0_+1]();
        dN_Poisson = new int [ct_+1]();

        double frac, frac0;
        double S_max = min(flux.flux_max, flux.S[flux.Num_Flux()-1]);

        frac = (log(S_max) - log(flux.S[ct_-1])) / (log(flux.S[ct_]) - log(flux.S[ct_-1]));
        if (ct0_ != 0) 
            frac0 = (log(flux.flux_min) - log(flux.S[ct0_-1])) / (log(flux.S[ct0_]) - log(flux.S[ct0_-1]));

        for (int i=0; i<ct_; ++i) dN[i] = flux.dNdlogS[i] * flux.dlogS * fourpi;
        dN[ct_] = flux.dNdlogS[ct_] * flux.dlogS * fourpi * frac;

        if (ct0_ == 1) {
            dN0[0] = flux.dNdlogS[0] * flux.dlogS * fourpi * frac0;
            dN[0] -= dN0[0]; }
        else if (ct0_ > 1) {
            for (int i=0; i<ct0_; ++i) dN0[i] = flux.dNdlogS[i] * flux.dlogS * fourpi;
            dN0[ct0_] = flux.dNdlogS[ct0_] * flux.dlogS * fourpi * frac0;

            for (int i=0; i<=ct_; ++i) dN[i] -= dN0[i];
        }

        if (map_dp == nullptr) map_dp = new double [npix];
    } else {
        cout << "dN has been calculated \n";
    }

    init = true;
}

void PoissonN::GenPN(rngHandle &rng) {

    viRngPoissonV(rng.vsl_method_poisson, rng.vsl_stream_poisson, ct_-ct0_+1, &dN_Poisson[ct0_], &dN[ct0_]);
    if (pixel_list == nullptr) {
        int num_uniform_max = *std::max_element(dN_Poisson, dN_Poisson+ct_);
        pixel_list = new int [2*num_uniform_max];
    }
}





template<typename T> void create_poisson_map(rngHandle &rng, FluxDensityInfo &flux, Healpix_Map<T> &map) {

    int ipix;

#pragma omp parallel shared (map) private (ipix) 
{
#pragma omp for
    for (ipix=0; ipix<map.Npix(); ++ipix) map[ipix] = 0.00e0;
}

    PoissonN poisson_number (rng, flux, map.Npix());

    double nu = 150.00E9;
    double Tcmb = 2.7255E0;
    double x = hPlanck * nu / (kBoltzmann * Tcmb);
    double dbnu = 2.00E0*kBoltzmann * (pow(nu,2.00E0)/pow(speedOfLight,2.00E0)) * (pow(x,2.00E0)*exp(x)) / pow((exp(x)-1.00E0), 2.00E0);
    double flx = 1.00E20 * dbnu; // 1000.0E0 * 1.00E26 * dbnu / 1.00E9

    double Omega_Pix = fourpi / ((double) map.Npix());
    double norm = flx * Omega_Pix;

    if (poisson_number.pixel_list == nullptr) {
        cout << "PoissonN initialization error \n";
        exit(1);
    }

    for (int i=0; i<poisson_number.Num_Element(); ++i) {
        if (poisson_number.dN_Poisson[i] == 0) continue;

        viRngUniform(rng.vsl_method_uniform, rng.vsl_stream_uniform, poisson_number.dN_Poisson[i], 
                     poisson_number.pixel_list, 0, map.Npix());

#pragma omp parallel shared (map) private (ipix)
{
#pragma omp for
        for (ipix=0; ipix<poisson_number.dN_Poisson[i]; ++ipix) 
            poisson_number.map_dp[ poisson_number.pixel_list[ipix] ] += flux.S[i];
}
    }

#pragma omp parallel shared (map) private (ipix)
{
#pragma omp for
    for (ipix=0; ipix<map.Npix(); ++ipix) map[ipix] = poisson_number.map_dp[ipix] / norm;
}
    

}

template void create_poisson_map(rngHandle &rng, FluxDensityInfo &flux, Healpix_Map<double> &map);
template void create_poisson_map(rngHandle &rng, FluxDensityInfo &flux, Healpix_Map<float> &map);


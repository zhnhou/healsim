#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

#include "healpix_map.h"
#include "healsim_rng.h"
#include "poisson_process.h"


FluxDensityInfo::FluxDensityInfo(string &FileName) {

    string      line;
    ifstream    infile(FileName.c_str());

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

template<typename T> void create_poission_map(rngHandle rng, FluxDensityInfo flux, Healpix_Map<T> map) {

}

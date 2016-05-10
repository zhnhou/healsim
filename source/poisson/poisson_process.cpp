#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "poisson_process.h"

using namespace std;

FluxDensity::FluxDensity(string &FileName) {

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
    Set(num_row);
}

FluxDensity::~FluxDensity(){}

void FluxDensity::assertArraySizes() const {
    planck_assert(multiequal( num_flux, S_.size(), S2p5dNdS_.size(), dNdS_.size(), 
    logS_.size(), dNdlogS_.size() ), "incorrect array sizes");

}

void FluxDensity::Set(int nflux) {
    num_flux = nflux;
    S_.alloc(num_flux);
    S2p5dNdS_.alloc(num_flux);
    dNdS_.alloc(num_flux);
    logS_.alloc(num_flux);
    dNdlogS_.alloc(num_flux);
}

//template<typename T> void create_poission_flux_density(const string FileName, P)

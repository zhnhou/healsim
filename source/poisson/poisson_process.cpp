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

}

FluxDensity::~FluxDensity(){}

//template<typename T> void create_poission_flux_density(const string FileName, P)

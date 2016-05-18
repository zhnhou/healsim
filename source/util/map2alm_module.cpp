#include <iostream>
#include <string>

#include "xcomplex.h"
#include "alm.h"
#include "healpix_map.h"
#include "alm_healpix_tools.h"
#include "healpix_data_io.h"
#include "map2alm_module.h"

using namespace std;

template<typename T> void map2alm_dummy (Healpix_Map<T> &map, Alm<xcomplex<T> > &alm, int lmax, int iter) {

    int nside = map.Nside();
    
    arr<double> weight;
    read_weight_ring ( string(getenv("HEALPIX")), nside, weight);

    alm(lmax, lmax);
    if (map.Scheme()==NEST) map.swap_scheme();

    map2alm_iter(map, alm, iter, weight);
}

template void map2alm_dummy (Healpix_Map<double> &map, Alm<xcomplex<double> > &alm, int lmax, int iter);
template void map2alm_dummy (Healpix_Map<float> &map, Alm<xcomplex<float> > &alm, int lmax, int iter);

#include "xcomplex.h"
#include "alm.h"
#include "healpix_map.h"
#include "healpix_data_io.h"

using namespace std;

template<typename T> int map2alm_dummy (Healpix_Map<T> &map, Alm<xcomplex<T> > &alm, int lmax, int iter=0) {

    int nside = map.Nside();
    
    arr<double> weight;
    get_weight_ring ( string(getenv("HEALPIX")), nside, weight);

    Alm<xcomplex<T> > alm(lmax, lmax);
    if (map.Scheme()==NEST) map.swap_scheme();

    map2alm_iter(map, alm, iter, weight);

    return 0;
}

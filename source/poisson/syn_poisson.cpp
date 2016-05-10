#include "error_handling.h"
#include "facilities.h"

int main (int argc, const char **argv) {
    PLANCK_DIAGNOSIS_BEGIN
    syn_poisson_module (argc, argv);
    PLANCK_DIAGNOSIS_END
}

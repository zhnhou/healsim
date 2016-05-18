CC = icpc

HPXCXX = $(HEALPIX)/src/cxx/linux_icc

CCFLAG = -O3 -openmp -std=c++11
LDFLAG = -O3 -openmp -std=c++11 -L$(CFITSIO)/lib -L$(HPXCXX)/lib -lhealpix_cxx -lcxxsupport -lsharp -lfftpack -lc_utils -lcfitsio -mkl

INCFLAG = -I./poisson/ -I./util/ -I$(CFITSIO)/include -I$(HPXCXX)/include -I$(MKLROOT)/include -I../include

OBJ_UTIL = util/map2alm_module.o

PRJ = syn_poisson
OBJ_POIS = poisson/healsim_rng.o poisson/poisson_process.o poisson/$(PRJ)_module.o

OBJECT = $(OBJ_UTIL) $(OBJ_POIS)


$(PRJ): $(PRJ).o $(OBJECT)
	$(CC) -o $(PRJ) $(PRJ).o $(OBJECT) $(LDFLAG) $(INCFLAG) -O3

$(PRJ).o: $(PRJ).cpp $(OBJECT)
	$(CC) $(CCFLAG) $(INCFLAG) -c $(PRJ).cpp

%.o: %.cpp
	$(CC) $(CCFLAG) $(INCFLAG) -c -o $*.o $*.cpp

clean:
	rm -rf *.o $(PRJ)

CC = icpc

HPXCXX = $(HEALPIX)/src/cxx/linux_icc

CCFLAG = -O3 -openmp -std=c++11
LDFLAG = -O3 -openmp -std=c++11 -L$(CFITSIO)/lib -L$(HPXCXX)/lib -lhealpix_cxx -lcxxsupport -lsharp -lfftpack -lc_utils -lcfitsio -mkl

INCFLAG = -I./ -I$(CFITSIO)/include -I$(HPXCXX)/include -I$(MKLROOT)/include -I../../include

PRJ = syn_poisson
OBJ = ../util/map2alm_module.o healsim_rng.o poisson_process.o $(PRJ)_module.o $(PRJ).o

%.o: %.cpp
	$(CC) $(CCFLAG) $(INCFLAG) -c $*.cpp

$(PRJ): $(OBJ)
	$(CC) -o $(PRJ) $(OBJ) $(LDFLAG) $(INCFLAG) -O3

clean:
	rm -rf *.o $(PRJ)

CC = icpc

HPXCXX = $(HEALPIX)/src/cxx/linux_icc

CCFLAG = -std=c++14 -O3 -openmp
LDFLAG = -std=c++14 -O3 -openmp -L$(HPXCXX)/lib -lcxxsupport -lhealpix_cxx -mkl

INCFLAG = -I./ -I$(HPXCXX)/include -I$(MKLROOT)/include -I../../include

PRJ = syn_poisson
OBJ = healsim_rng.o poisson_process.o $(PRJ)_module.o $(PRJ).o

%.o: %.cpp
	$(CC) $(CCFLAG) $(INCFLAG) -c $*.cpp

$(PRJ): $(OBJ)
	$(CC) -o $(PRJ) $(OBJ) $(LDFLAG) $(INCFLAG) -O3

clean:
	rm -rf *.o $(PRJ)

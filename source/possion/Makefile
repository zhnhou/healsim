CC = icpc

HPXCXX = $(HEALPIX)/src/cxx/linux_icc

CCFLAG = -O3 -openmp
LDFLAG = -O3 -openmp -L$(HPXCXX)/lib -lcxxsupport -lhealpix_cxx

INCFLAG = -I$(HPXCXX)/include -I../../include

PRJ = syn_poisson
OBJ = $(PRJ)_module.o $(PRJ).o

%.o: %.cpp
	$(CC) $(CCFLAG) $(INCFLAG) -c $*.cpp

$(PRJ): $(OBJ)
	$(CC) -o $(PRJ) $(OBJ) $(LDFLAG) $(INCFLAG) -O3

clean:
	rm -rf *.o $(PRJ)

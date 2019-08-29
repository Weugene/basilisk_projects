QCC_COMPILER = $(BASILISK)/qcc
CC99="/opt/intel/compilers_and_libraries_2018.2.199/linux/mpi/intel64/bin//mpicc -std=c99"
CFLAGS=-DDUMP=1 -Wall -O2
CLIBS=-lm -L$(BASILISK)/gl  $(OPENGLIBS)
CINCL=-I$(BASILISK)

alphaT:
	CC99='mpicc -std=c99' $(QCC_COMPILER) $(CFLAGS)  -D_MPI=1  alphaT.c $(CINCL) $(CLIBS)

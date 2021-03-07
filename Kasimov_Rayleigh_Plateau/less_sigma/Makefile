#CFLAGS += -O2
#include $(BASILISK)/Makefile.defs
QCC_COMPILER = qcc
CC99='/opt/intel/compilers_and_libraries_2018.2.199/linux/mpi/intel64/bin//mpicc -std=c99'
Kasimov_RP:
#	$(QCC_COMPILER) -Wall -O2 Kasimov.c -L$(BASILISK)/gl $(OPENGLIBS) -lm
#	./a.out > out 2> log
#	animate out.ppm
	CC99='mpicc -std=c99' $(QCC_COMPILER)   -O2 -D_MPI=1  Kasimov.c  -L$(BASILISK)/gl $(OPENGLIBS) -lm
	mpirun -np 10 ./a.out > out 2> log
#	animate movie.mp4.ppm
#DO NOT GENERATE VTK IN PARALLEL
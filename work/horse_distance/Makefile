#CFLAGS += -O2
#include $(BASILISK)/Makefile.defs
QCC_COMPILER = qcc
CC99="/opt/intel/compilers_and_libraries_2018.2.199/linux/mpi/intel64/bin//mpicc -std=c99"
distance:
	$(QCC_COMPILER) -Wall -O2 distance.c -lm  -L$(BASILISK)/gl $(OPENGLIBS)
	./a.out > out 2> log
#	animate out.ppm
parallel_distance:
	CC99='mpicc -std=c99' $(QCC_COMPILER) -Wall  -O2 -D_MPI=1  distance.c -lm -L$(BASILISK)/gl $(OPENGLIBS)
	mpirun -np 10 ./a.out > out_ 2> log_
#	animate out.ppm

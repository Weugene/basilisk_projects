#CFLAGS += -O2
#include $(BASILISK)/Makefile.defs
QCC_COMPILER = /home/weugene/basilisk/src/qcc
CC99='/opt/intel/compilers_and_libraries_2018.2.199/linux/mpi/intel64/bin//mpicc -std=c99'

brusselator:
#	$(QCC_COMPILER) -Wall -O2 brusselator.c -lm
#	./a.out > out 2> log
#	animate out.ppm
	CC99='mpicc -std=c99' $(QCC_COMPILER)   -O2 -D_MPI=1  brusselator.c -lm -L$(BASILISK)/gl $(OPENGLIB)
	mpirun -np 8 ./a.out > out 2> log
#	animate out
brusselator3samples:
#	$(QCC_COMPILER) -Wall -O2 brusselator.c -lm
#	./a.out > out 2> log
#	animate out.ppm
	CC99='mpicc -std=c99' $(QCC_COMPILER)   -O2 -D_MPI=1  brusselator.c -lm -L$(BASILISK)/gl $(OPENGLIB)
	mpirun -np 8 ./a.out > out_ 2> log_
#	animate out
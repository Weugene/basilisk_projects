#CFLAGS += -O2
#include $(BASILISK)/Makefile.defs
#BASILISK=/home/weugene/basilisk/src
#QCC_COMPILER = $(BASILISK)/qcc
#CC99="/opt/intel/compilers_and_libraries_2018.2.199/linux/mpi/intel64/bin//mpicc -std=c99"
CFLAGS=-DDUMP=1 -Wall -O2
CLIBS=-lm -L$(BASILISK)/gl  $(OPENGLIB)
CINCL=-I$(BASILISK)
CALC_TIME=72000

jet_dischargejose:
	qcc   $(CFLAGS) jet_dischargejose.c $(CINCL) $(CLIBS)
	./a.out > out 2> log
#	mvlc movie.mp4
#	CC99='mpicc -std=c99' qcc $(CFLAGS)  -D_MPI=1  jet_dischargejose.c $(CINCL) $(CLIBS)
#	mpirun -np 5 ./a.out
#	mvlc movie.mp4


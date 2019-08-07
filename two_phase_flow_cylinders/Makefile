#CFLAGS += -O2
#include $(BASILISK)/Makefile.defs
QCC_COMPILER = /home/weugene/basilisk/src/qcc
CC99="/opt/intel/compilers_and_libraries_2018.2.199/linux/mpi/intel64/bin//mpicc -std=c99"
CFLAGS=-DDUMP=1 -Wall -O2
CLIBS=-lm -L$(BASILISK)/gl  $(OPENGLIBS)
CINCL=-I$(BASILISK)
CALC_TIME=7200000

two_phase_flow_cylinders:
#	qcc   $(CFLAGS) two_phase_flow_cylinders.c $(CINCL) $(CLIBS)
#	./a.out > out 2> log
	CC99='mpicc -std=c99' qcc $(CFLAGS)  -D_MPI=1  two_phase_flow_cylinders.c $(CINCL) $(CLIBS)
#	mpirun -np 20 ./a.out -m $(CALC_TIME) 12 0.4 > out 2> log

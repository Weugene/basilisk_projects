#CFLAGS += -O2
#include $(BASILISK)/Makefile.defs
BASILISK=/home/weugene/basilisk/src
QCC_COMPILER = $(BASILISK)/qcc
CC99="/opt/intel/compilers_and_libraries_2018.2.199/linux/mpi/intel64/bin//mpicc -std=c99"
CFLAGS=-DDUMP=1 -Wall -O2
CLIBS=-lm -L$(BASILISK)/gl  $(OPENGLIB)
CINCL=-I$(BASILISK)
CALC_TIME=36000

mask:
	qcc   $(CFLAGS) mask.c $(CINCL) $(CLIBS)
	./a.out


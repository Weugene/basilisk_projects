#CFLAGS += -O2
#include $(BASILISK)/Makefile.defs
QCC_COMPILER = /home/weugene/basilisk/src/qcc
CC99='/opt/intel/compilers_and_libraries_2018.2.199/linux/mpi/intel64/bin//mpicc -std=c99'
all:
	CC99='mpicc -std=c99' $(QCC_COMPILER) -Wall -O2 -D_MPI=1  explosion.c -lm && mpirun -np 8 ./a.out > out 2> log && animate out
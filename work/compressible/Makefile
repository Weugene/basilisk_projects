#CFLAGS += -O2
#include $(BASILISK)/Makefile.defs
QCC_COMPILER = /home/weugene/basilisk/src/qcc
CC99='/opt/intel/compilers_and_libraries_2018.2.199/linux/mpi/intel64/bin//mpicc -std=c99'
all:
	$(QCC_COMPILER) -Wall -O2 explosion.c -lm
	./a.out > out.ppm 2> log
	animate out.ppm
	#CC99='mpicc -std=c11' $(QCC_COMPILER) -Wall -O2 -D_MPI=1  explosion.c -lm && mpirun -np 8 ./a.out > out_ 2> log && animate out
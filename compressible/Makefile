#CFLAGS += -O2
#include $(BASILISK)/Makefile.defs
QCC_COMPILER = qcc
CC99='/opt/intel/compilers_and_libraries_2018.2.199/linux/mpi/intel64/bin//mpicc -std=c99'
all:
	$(QCC_COMPILER) -Wall -O2 explosion.c -lm -L$(BASILISK)/gl
	./a.out > out 2> log
#	animate out.ppm
#	CC99='mpicc -std=c99' $(QCC_COMPILER) -Wall -O2 -D_MPI=1  explosion.c -L$(BASILISK)/gl $(OPENGLIBS)
#	mpirun -np 20 ./a.out > out_ 2> log
#	animate out.ppm
anim:
	animate out.ppm
anim_:
	animate out_

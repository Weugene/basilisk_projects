#CFLAGS += -O2
#include $(BASILISK)/Makefile.defs
QCC_COMPILER = /home/weugene/basilisk/src/qcc
CC99="/opt/intel/compilers_and_libraries_2018.2.199/linux/mpi/intel64/bin//mpicc -std=c99"
OPENGLIB="-lfb_glx -lGLU -lGLEW -lGL -lX11 -lfb_osmesa -lglutils"
flow_over_cylinder:
	$(QCC_COMPILER) -Wall -O2 flow_over_cylinder.c -lm -L$BASILISK/gl  $(OPENGLIB)
	./a.out > out 2> log
#	animate out.ppm
#	CC99='mpicc -std=c99' $(QCC_COMPILER) -O2 -D_MPI=1 flow_over_cylinder.c -lm -L$(BASILISK)/gl $(OPENGLIB)
#	mpirun -np 1 ./a.out > out 2> log
	animate vort.ppm && animate p.ppm
MPI_flow_over_cylinder:
	CC99='mpicc -std=c99' $(QCC_COMPILER) -Wall -D_MPI=1 -O2 flow_over_cylinder.c -lm -L$BASILISK/gl  $(OPENGLIB)
	mpirun -np 10 ./a.out > out_ 2> log_
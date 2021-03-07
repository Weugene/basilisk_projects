#CFLAGS += -O2
#include $(BASILISK)/Makefile.defs
QCC_COMPILER = /home/weugene/basilisk/src/qcc
CC99='/opt/intel/compilers_and_libraries_2018.2.199/linux/mpi/intel64/bin//mpicc -std=c99'
#-L$BASILISK/gl  -lGLU -lGLEW -lGL -lX11 #-lglutils
all:
	$(QCC_COMPILER)  -Wall -O2    curvature.c -lm -lvofi
	./a.out > out 2> log
#	animate f.ppm
#	animate grid.ppm
#	CC99='mpicc -std=c99' $(QCC_COMPILER) -Wall -O2 -D_MPI=1 -L$BASILISK/gl -lglutils -lfb_glx -lGLU -lGLEW -lGL -lX11 explosion.c
#	mpirun -np 20 ./a.out > out_ 2> log
#	animate out.ppm
anim:
	animate f.ppm
	animate grid.ppm
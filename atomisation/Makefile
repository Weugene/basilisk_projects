#CFLAGS += -O2
#include $(BASILISK)/Makefile.defs
#QCC_COMPILER = /home/weugene/basilisk/src/qcc
#CC99='/opt/intel/compilers_and_libraries_2018.2.199/linux/mpi/intel64/bin//mpicc -std=c99'
atomisation:
	qcc -Wall -O2 atomisation.c -L$(BASILISK)/gl -L/opt/local/lib -lm -lfb_osmesa -lOSMesa -lglutils -framework OpenGL -framework GLUT
	./a.out > out.ppm 2> log
#	animate out.ppm
#	CC99='mpicc -std=c99' qcc   -O2 -D_MPI=1  atomisation.c  -L$(BASILISK)/gl $(OPENGLIBS) -lm
#	mpirun -np 10 ./a.out > out 2> log
#	animate movie.mp4.ppm

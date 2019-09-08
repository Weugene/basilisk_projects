CFLAGS += -DDUMP=1 -Wall -O2 -g
CLIBS=-lm -L$(BASILISK)/gl  -L/opt/local/lib  -lfb_osmesa -lOSMesa -lglutils -framework OpenGL -framework GLUT
#-lfb_osmesa -lglutils
#$(OPENGLIBS) -lX11 -lGLU

karman:
#	qcc   $(CFLAGS) karman.c $(CINCL) $(CLIBS)
#	./a.out > out 2> log
#	mvlc movie.mp4
	CC99='mpicc -std=c99' qcc $(CFLAGS)  -D_MPI=1  karman.c  $(CLIBS)
	mpirun -np 4 ./a.out  > out 2> log
#	mvlc movie.mp4

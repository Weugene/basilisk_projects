CFLAGS=-DDUMP=1 -Wall -O2 -events
CLIBS=-lm $(OPENGLIBS)
CINCL=-I$(BASILISK)
CORES=1
CFILE=bone
bone3D:
	qcc   $(CFLAGS) $(CFILE).c $(CINCL) $(CLIBS)
#	./a.out
#	CC99='mpicc -std=c99' qcc $(CFLAGS)  -D_MPI=1  $(CFILE).c $(CINCL) $(CLIBS)
#	mpirun -np $(CORES) ./a.out > out 2> log
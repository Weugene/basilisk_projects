CFLAGS=-DDUMP=1 -Wall -O2 -events
CLIBS=-lm $(OPENGLIBS)
CINCL=-I$(BASILISK)
CFILE=2phase_BP
2phase_BP:
#	qcc $(CFLAGS) $(CFILE).c $(CINCL) $(CLIBS)
#	./a.out
#> out 2> log
	CC99='mpicc -std=c99' qcc $(CFLAGS)  -D_MPI=1  $(CFILE).c $(CINCL) $(CLIBS)
	mpirun -np 4 ./a.out  > out 2> log
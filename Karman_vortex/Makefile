CFLAGS=-DDUMP=1 -Wall -O2 -events
CLIBS=-lm $(OPENGLIBS)
CFILE=karman_vortex
karman_vortex:
	rm -rf .qcc* kk* a.out dump* out log
	$(BASILISK)/qcc   $(CFLAGS) $(CFILE).c $(CINCL) $(CLIBS)
	./a.out
#	CC99='mpicc -std=c99' qcc $(CFLAGS)  -D_MPI=1  $(CFILE).c  $(CLIBS)
#	mpirun -np 4 ./a.out  > out 2> log


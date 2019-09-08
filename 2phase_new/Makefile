CFLAGS=-DDUMP=1 -Wall -O2
CLIBS=-lm $(OPENGLIBS)
CINCL=-I$(BASILISK)

two_phase_new:
#	qcc   $(CFLAGS) two_phase_new.c $(CINCL) $(CLIBS)
#	./a.out > out 2> log
#	mvlc movie.mp4
	CC99='mpicc -std=c99' qcc $(CFLAGS)  -D_MPI=1  two_phase_new.c $(CINCL) $(CLIBS)
#	mpirun -np 10 ./a.out  > out 2> log
#	mvlc movie.mp4
#two_phase_new.c
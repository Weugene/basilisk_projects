CFLAGS=-DDUMP=1 -Wall -O2
CLIBS=-lm -L$(BASILISK)/gl  $(OPENGLIBS)
CALC_TIME=72000

tangaroa:
#	qcc   $(CFLAGS) tangaroa.c $(CINCL) $(CLIBS)
#	./a.out > out 2> log
#	mvlc movie.mp4
	CC99='mpicc -std=c99' qcc $(CFLAGS)  -D_MPI=1  tangaroa.c $(CINCL) $(CLIBS)
	mpirun -np 4 ./a.out -m $(CALC_TIME) 9 0.4 > out 2> log
	mvlc movie.mp4


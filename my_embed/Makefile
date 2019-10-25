CFLAGS=-DDUMP=1 -Wall -O2 -events
CLIBS=-lm $(OPENGLIBS)
CFILE=my_embed
my_embed:
	rm -rf .qcc* br* a.out dump* out log
	$(BASILISK)/qcc   $(CFLAGS) $(CFILE).c $(CINCL) $(CLIBS)
	./a.out
#	 > out 2> log &
#	mvlc movie.mp4
#	CC99='mpicc -std=c99' qcc $(CFLAGS)  -D_MPI=1  $(CFILE).c  $(CLIBS)
#	mpirun -np 4 ./a.out  > out 2> log
#	mvlc movie.mp4
#   qcc   $(CFLAGS) $(CFILE).c $(CINCL) $(CLIBS)
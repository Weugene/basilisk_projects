CFLAGS=-DDUMP=1 -Wall -O2
CLIBS=-lm $(OPENGLIBS)
CFILE=wannier_bk
wannier_bk:
	rm -rf .qcc* wannier_* a.out dump* out log
	$(BASILISK)/qcc   $(CFLAGS) $(CFILE).c $(CINCL) $(CLIBS)
	./a.out > out 2> log
#	CC99='mpicc -std=c99' qcc $(CFLAGS)  -D_MPI=1  $(CFILE).c  $(CLIBS)
#	mpirun -np 4 ./a.out  > out 2> log

#   qcc   $(CFLAGS) $(CFILE).c $(CINCL) $(CLIBS)

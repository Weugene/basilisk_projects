CFILE=hjump
hjump:
	rm -rf log out hjump_* snapshot-* .qcc*
	qcc -Wall -O2 $(CFILE).c -lm $(OPENGLIBS)
	./a.out > out 2> log
#	CC99='mpicc -std=c99' qcc  -Wunused-result -O2 -D_MPI=1  $(CFILE).c -lm  $(OPENGLIBS)
#	mpirun -np 4 ./a.out > out 2> log

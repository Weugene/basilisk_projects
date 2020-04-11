CFILE=vof_gi
vof_gi:
	rm -rf .qcc* log out rk_* snapshot-* nameYouWant*
	qcc -O2 -events $(CFILE).c -lm $(OPENGLIBS)
	./a.out 10 200 1 1 1 10> out 2> log
#	CC99='mpicc -std=c99' qcc -events -Wunused-result -O2 -D_MPI=1  $(CFILE).c -lm $(OPENGLIBS)
#	mpirun -np 4 ./a.out > out 2> log

#	qcc -O2 -events vof_gi.c -lm $OPENGLIBS
#	CC99='mpicc -std=c99' qcc -events -Wunused-result -O2 -D_MPI=1  vof_gi.c -lm $OPENGLIBS

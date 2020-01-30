CFILE=roughness
roughness:
	rm -rf log out porous_BP_*
	qcc -Wall -O2 $(CFILE).c -lm -events $(OPENGLIBS)
	./a.out > out 2> log
#	CC99='mpicc -std=c99' qcc  -Wunused-result -O2 -D_MPI=1  $(CFILE).c -lm  $(OPENGLIBS)
#	mpirun -np 10 ./a.out > out 2> log

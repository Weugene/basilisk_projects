CFILE=atomisation_BP
atomisation_BP:
	rm -rf .qcc* hr_* atom_* dump out* log* a.out
#	qcc -Wall -O2 $(CFILE).c -lm -events -L$(BASILISK)/gl  $(OPENGLIBS)
#	./a.out >out 2>log
	CC99='mpicc -std=c99' qcc  -Wunused-result -O2 -D_MPI=1  $(CFILE).c -lm -L$(BASILISK)/gl $(OPENGLIBS)
	mpirun -np 4 ./a.out

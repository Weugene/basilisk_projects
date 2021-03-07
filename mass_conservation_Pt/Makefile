CFILE=mass_cons_Pt
mass_cons_Pt:
	rm -rf .qcc* log out mc_* snapshot-* nameYouWant*
	qcc -Wall -O2 -events $(CFILE).c -lm $(OPENGLIBS)
	./a.out > out 2> log
#	CC99='mpicc -std=c99' qcc -Wunused-result -O2 -D_MPI=1 -events $(CFILE).c -lm $(OPENGLIBS)
#	mpirun -np 4 ./a.out > out 2> log

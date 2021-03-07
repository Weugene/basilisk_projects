porous:
	qcc -Wall -O2 porous.c -lm -L$(BASILISK)/gl  $(OPENGLIBS)
	./a.out > out 2> log
#	animate out.ppm
#	CC99='mpicc -std=c99' qcc  -Wunused-result -O2 -D_MPI=1  porous.c -lm -L$(BASILISK)/gl $(OPENGLIB)
#	mpirun -np 10 ./a.out > out 2> log
#	animate out.ppm

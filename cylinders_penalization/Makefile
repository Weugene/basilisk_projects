
cylinders_penalization:
	rm -rf .qcc* hr_* br_* dump out* log* a.out
#	qcc -Wall -O2 cylinders_penalization.c -lm -events -L$(BASILISK)/gl  $(OPENGLIBS)
#	./a.out >out 2>log
	CC99='mpicc -std=c99' qcc  -Wunused-result -O2 -D_MPI=1  cylinders_penalization.c -lm -L$(BASILISK)/gl $(OPENGLIBS)
	mpirun -np 4 ./a.out

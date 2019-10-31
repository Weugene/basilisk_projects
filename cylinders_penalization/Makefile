cylinders_penalization:
	rm -rf .qcc* hr_* br_* dump out* log* a.out
	qcc -Wall -O2 cylinders_penalization.c -lm -L$(BASILISK)/gl  $(OPENGLIBS)
	./a.out	> out 2> log
#	animate out.ppm
#	CC99='mpicc -std=c99' $(QCC_COMPILER)  -Wunused-result -O2 -D_MPI=1  cylinders_test.c -lm -L$(BASILISK)/gl $(OPENGLIB)
#	mpirun -np 10 ./a.out > out 2> log
#	animate out.ppm

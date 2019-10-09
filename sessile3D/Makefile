
DIMENSION=""
#DIMENSION="3D"
sessile:
	qcc -Wall -events -O2 sessile$(DIMENSION).c -lm $(OPENGLIB)
	./a.out > out 2> log
#	animate out.ppm
#	CC99='mpicc -std=c99' $(QCC_COMPILER) -Wall -O2 -D_MPI=1  sessile$(DIMENSION).c -lm $(OPENGLIBS)
#	mpirun -np 10 ./a.out > out 2> log
#	animate out.ppm

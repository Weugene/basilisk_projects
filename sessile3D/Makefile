
DIMENSION=3D
#DIMENSION="3D"
sessile:
	rm -rf .qcc* log out out_ ses_*
#	qcc -Wall -events -O2 sessile$(DIMENSION).c -lm $(OPENGLIBS)
#	./a.out > out 2> log
	CC99='mpicc -std=c99' qcc -Wall -O2 -D_MPI=1  sessile$(DIMENSION).c -lm $(OPENGLIBS)
	mpirun -np 5 ./a.out > out 2> log

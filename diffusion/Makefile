QCC_COMPILER = qcc
diffusion:
	$(QCC_COMPILER)  -O2    diffusion.c -lm
	./a.out > out 2> log
#	animate f.ppm
#	animate grid.ppm
#	CC99='mpicc -std=c99' $(QCC_COMPILER) -Wall -O2 -D_MPI=1 diffusion.c -lm -L$(BASILISK)/gl $(OPENGLIB)
#	mpirun -np 10 ./a.out > out_ 2> log
#	animate out.ppm
anim:
	animate f.ppm
	animate grid.ppm

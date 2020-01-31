CFILE=roughness
all: roughnessP

roughnessS:
	rm -rf log out rg_*
	qcc -Wall -O2 $(CFILE).c -lm -events $(OPENGLIBS)
	./a.out > out 2> log

roughnessP:
	rm -rf log out rg_*
	CC99='mpicc -std=c99' qcc  -Wunused-result -O2 -D_MPI=1  $(CFILE).c -lm  $(OPENGLIBS)
	mpirun -np 4 ./a.out > out 2> log

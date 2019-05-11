QCC_COMPILER = $(BASILISK)/qcc
#CC99="/opt/intel/compilers_and_libraries_2018.2.199/linux/mpi/intel64/bin//mpicc -std=c99"


three_phase:
	$(QCC_COMPILER) -Wall -O2 three_phase.c -lm -L$(BASILISK)/gl  $(OPENGLIBS)
	./a.out > out 2> log

#	CC99='mpicc -std=c99' $(QCC_COMPILER) -Wall -O2 -D_MPI=1  three_phase.c -lm -L$(BASILISK)/gl $(OPENGLIBS)
#	mpirun -np 10 ./a.out > out 2> log


QCC_COMPILER = $(BASILISK)/qcc
#CC99="/opt/intel/compilers_and_libraries_2018.2.199/linux/mpi/intel64/bin//mpicc -std=c99"


reversed_3f-implicit:
	$(QCC_COMPILER) -Wall -O2 reversed_3f-implicit.c -lm -L$(BASILISK)/gl  $(OPENGLIBS)
	./a.out > out 2> log

#	CC99='mpicc -std=c99' $(QCC_COMPILER) -Wall -O2 -D_MPI=1  reversed_3f-implicit.c -lm -L$(BASILISK)/gl $(OPENGLIBS)
#	mpirun -np 10 ./a.out > out 2> log


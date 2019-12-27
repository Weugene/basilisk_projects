tangaroa_BP:
#	rm -rf .qcc* hr_* br_* dump out* log* a.out
#	qcc -Wall -O2 -DDUMP=1 tangaroa_BP.c -lm  -L$(BASILISK)/gl  $(OPENGLIBS)
#	./a.out
# >out 2>log
	CC99='mpicc -std=c99' qcc -O2 -D_MPI=1 -DDUMP=1 -events tangaroa_BP.c -lm -L$(BASILISK)/gl $(OPENGLIBS)
	mpirun -np 10 ./a.out >out 2>log

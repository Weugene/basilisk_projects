CFILE=3D_composite
3D_composite:
	rm -rf .qcc* hr_* br_* dump out* log* a.out
	qcc -Wall -O2 -DDUMP=1 $(CFILE).c -lm  $(OPENGLIBS)
	./a.out >out 2>log
#	CC99='mpicc -std=c99' qcc -O2 -D_MPI=1 -DDUMP=1 $(CFILE).c -lm $(OPENGLIBS)
#	mpirun -np 10 ./a.out >out 2>log &

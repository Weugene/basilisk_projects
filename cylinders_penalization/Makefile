
cylinders_penalization:
	rm -rf .qcc* hr_* br_* dump out log a.out cyl_br* res/cyl_br*
	qcc -Wall -O2 -events cylinders_penalization.c -lm $(OPENGLIBS)
	./a.out > out_bp 2> log_bp
#	CC99='mpicc -std=c99' qcc  -Wunused-result -O2 -D_MPI=1  cylinders_penalization.c -lm -L$(BASILISK)/gl $(OPENGLIBS)
#	mpirun -np 4 ./a.out

CFLAG= -O2 -events
CFILE=tube_triangle.c
tube_triangle:
	rm -rf .qcc* log_bp out_bp snapshot-*
#	qcc  $(CFLAG) $(CFILE) -lm
#	./a.out > out_bp 2> log_bp
	CC99='mpicc -std=c99' qcc -O2 -D_MPI=1 -events $(CFILE) -lm
	mpirun -np 2 ./a.out  > out_bp 2> log_bp


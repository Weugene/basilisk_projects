CFLAG= -Wall -O2 -events -g
all:
	echo "choose tube_eb or tube_bp"
tube_eb:
	rm -rf .qcc* log out res/tube_eb_* snapshot-* tube_eb.pvd
	qcc  $(CFLAG) tube_eb.c -lm $(OPENGLIBS)
	./a.out > out_eb 2> log_eb
#	CC99='mpicc -std=c99' $(QCC_COMPILER)  -Wunused-result -O2 -D_MPI=1  tube_eb.c -lm $(OPENGLIBS)
#	mpirun -np 4 ./a.out > out_eb 2> log_eb
tube_bp:
	rm -rf .qcc* log_bp out_bp res/tube_bp_* snapshot-* tube_bp.pvd
#	qcc  $(CFLAG) tube_bp.c -lm $(OPENGLIBS)
#	./a.out > out_bp 2> log_bp
	CC99='mpicc -std=c99' qcc  -Wunused-result -O2 -D_MPI=1 -events tube_bp.c -lm $(OPENGLIBS)
	mpirun -np 4 ./a.out > out_bp 2> log_bp
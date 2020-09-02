CFLAG= -O2 -events
all:
	echo "choose tube_eb or tube_bp"
tube_eb:
	rm -rf .qcc* log out res/tube_eb_* snapshot-* tube_eb.pvd
	qcc  $(CFLAG) tube_eb.c -lm #$(OPENGLIBS)
	./a.out > out_eb 2> log_eb
#	CC99='mpicc -std=c99' $(QCC_COMPILER)  -Wunused-result -O2 -D_MPI=1  tube_eb.c -lm #$(OPENGLIBS)
#	mpirun -np 4 ./a.out > out_eb 2> log_eb
tube_bp:
	rm -rf .qcc* log_bp out_bp snapshot-*
#res/tube_bp_* tube_bp.pvd
#	qcc  $(CFLAG) tube_bp.c -lm
#$(OPENGLIBS)
#	./a.out > out_bp 2> log_bp
	CC99='mpicc -std=c99' qcc -O2 -D_MPI=1 -events tube_bp.c -lm
#$(OPENGLIBS)
#mpirun -np 4 ./a.out  > out_bp 2> log_bp
convertD2P_tar:
	rm -rf convertD2P
	CC99='mpicc -std=c99' qcc $(CFLAG) -D_MPI=1 convertD2P.c -lm -o convertD2P
	mpirun -np 4 ./convertD2P 10 30 1 0
#	qcc  $(CFLAG) convertD2P.c -lm -o convertD2P
#	./convertD2P 10 30 1 0
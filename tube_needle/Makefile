CFLAG= -O2 -events
all:
	echo 'tube_bp'
tube_bp:
	rm -rf .qcc* log_bp out_bp snapshot-*
	CC99='mpicc -std=c99' qcc -O2 -D_MPI=1 -events tube_bp.c -lm
#	mpirun -np 2 ./a.out
convertD2P_tar:
	rm -rf convertD2P
	CC99='mpicc -std=c99' qcc $(CFLAG) -D_MPI=1 convertD2P.c -lm -o convertD2P
	mpirun -np 4 ./convertD2P 10 30 1 0
#	qcc  $(CFLAG) convertD2P.c -lm -o convertD2P
#	./convertD2P 10 30 1 0
convert_single_tar:
	rm -rf convert_single
#    CC99='mpicc -std=c99' qcc $(CFLAG) -D_MPI=1 convertD2P_single.c -lm -o convert_single_p
#    mpirun -np 4 ./convert_single_p 10 30 1 0
	qcc  $(CFLAG) convertD2P_single.c -lm -o convert_single
#    ./convert_single 10 30 1 0

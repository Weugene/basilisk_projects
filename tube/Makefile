CFLAG= -O2 -events
ACTIVE_ENV=source ~/.bash_profile

all:
	echo "choose tube_eb or tube_bp"
tube_eb:
	rm -rf .qcc* log out res/tube_eb_* snapshot-* tube_eb.pvd
	qcc  $(CFLAG) tube_eb.c -lm #$(OPENGLIBS)
#	./a.out > out_eb 2> log_eb
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
#	mpirun -np 4 ./convertD2P 10 30 1 0
#	qcc  $(CFLAG) convertD2P.c -lm -o convertD2P
#	./convertD2P 10 30 1 0
convert_single_tar:
	rm -rf convert_single
#    CC99='mpicc -std=c99' qcc $(CFLAG) -D_MPI=1 convertD2P_single.c -lm -o convert_single_p
#    mpirun -np 4 ./convert_single_p 10 30 1 0
	qcc  $(CFLAG) convertD2P_single.c -lm -o convert_single
#    ./convert_single 10 30 1 0
convert_reduce_maxl:
	rm -rf convert_reduce_maxlevel
	qcc  $(CFLAG) convert_single_reduce_maxlevel.c -lm -g -o convert_reduce_maxlevel
#    ./convert_reduce_maxlevel 10 30 1 0
tube_bp_rel:
#	qcc  $(CFLAG) tube_bp_rel.c -lm
#	./a.out > out_bp 2> log_bp
	CC99='mpicc -std=c99' qcc -O2 -D_MPI=1 -events tube_bp_rel.c -lm
#	mpirun -np 8 ./a.out
#   > out_bp 2> log_bp
tube_bp_rel_buggy:
#	qcc  $(CFLAG) tube_bp_rel_buggy.c -lm
#	./a.out > out_bp 2> log_bp
	CC99='mpicc -std=c99' qcc -O2 -D_MPI=1 -events tube_bp_rel_buggy.c -lm
#	mpirun -np 8 ./a.out
#   > out_bp 2> log_bp
tube_bp_rel_server:
	path="$(HOME)/wbasilisk/tube"
	rsync -av tube_bp_rel.c Makefile h:$(path)
#	ssh h "source /home/e.sharaborin/.bash_profile && cd  ~/wbasilisk/tube/ && make tube_bp_rel &"
#	ssh h "cd ~/wbasilisk/tube/ && tail -n 20 log_bp"
#	clear && echo -en '\''\e[3J'\'''
#	ssh h "while sleep 5; do cd ~/wbasilisk/tube/; tail -n 5 log_bp; done &"

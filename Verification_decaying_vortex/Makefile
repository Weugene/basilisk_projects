CFLAG=-Wno-logical-op-parentheses -Wall -O2 -events -g
#-catch -g
all:
	echo "choose decaying_vortex_EB or decaying_vortex_Popinet or decaying_vortex_BP"
decaying_vortex_EB:
	rm -rf .qcc* log out res/vortex_EB_* snapshot-*
	qcc  $(CFLAG) decaying_vortex_EB.c -lm $(OPENGLIBS)
	./a.out > out_EB 2> log_EB
#	CC99='mpicc -std=c99' qcc -events -O2 -D_MPI=1  $(CFILE).c -lm $(OPENGLIBS)
#	mpirun -np 4 ./a.out > out 2> log
decaying_vortex_Popinet:
	rm -rf .qcc* log out res/vortex_Popinet_* snapshot-* nameYouWant*
	qcc $(CFLAG)  decaying_vortex_Popinet.c -lm $(OPENGLIBS)
	./a.out > out_P 2> log_P
decaying_vortex_BP:
	rm -rf .qcc* log out snapshot-* nameYouWant*
	qcc $(CFLAG) decaying_vortex_BP.c -lm $(OPENGLIBS)
	./a.out > out_BP 2> log_BP
#	CC99='mpicc -std=c99' qcc -Wall -O2 -DDUMP=1 -D_MPI=1 -DDIVRHS=0 -events decaying_vortex_BP.c -lm $(OPENGLIBS) -o parallel.out
#	mpirun -np 4 ./parallel.out > out_BP 2> log_BP
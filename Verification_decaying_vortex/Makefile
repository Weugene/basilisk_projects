all:
	echo "choose decaying_vortex_EB or starting"
decaying_vortex_EB:
	rm -rf .qcc* log out res/vortex_EB_* snapshot-*
	qcc -Wall -O2 -events -Wunused-result -Wlogical-op-parentheses decaying_vortex_EB.c -lm $(OPENGLIBS)
	./a.out > out_EB 2> log_EB
#	CC99='mpicc -std=c99' qcc -events -Wunused-result -Wlogical-op-parentheses -O2 -D_MPI=1  $(CFILE).c -lm $(OPENGLIBS)
#	mpirun -np 4 ./a.out > out 2> log
decaying_vortex_Popinet:
	rm -rf .qcc* log out res/vortex_Popinet_* snapshot-* nameYouWant*
	qcc -Wall -O2 -g -events -Wunused-result -Wlogical-op-parentheses decaying_vortex_Popinet.c -lm $(OPENGLIBS)
	./a.out > out_P 2> log_P
decaying_vortex_BP:
	rm -rf .qcc* log out res/vortex_BP_* snapshot-* nameYouWant*
	qcc -Wall -O2 -events -Wunused-result -Wlogical-op-parentheses decaying_vortex_BP.c -lm $(OPENGLIBS)
	./a.out > out_BP 2> log_BP
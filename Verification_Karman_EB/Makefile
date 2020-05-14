
all:
	echo "choose karman_EB or starting"
karman_EB:
	rm -rf .qcc* log out rkEB_* snapshot-* nameYouWant*
	qcc -Wall -O2 -events -Wunused-result -Wlogical-op-parentheses karman_EB.c -lm $(OPENGLIBS)
	./a.out > out 2> log
#	CC99='mpicc -std=c99' qcc -events -Wunused-result -Wlogical-op-parentheses -O2 -D_MPI=1  $(CFILE).c -lm $(OPENGLIBS)
#	mpirun -np 4 ./a.out > out 2> log
starting:
	rm -rf .qcc* log out start_* snapshot-* nameYouWant*
	qcc -Wall -O2 -events -Wunused-result -Wlogical-op-parentheses starting.c -lm $(OPENGLIBS)
	./a.out > out 2> log
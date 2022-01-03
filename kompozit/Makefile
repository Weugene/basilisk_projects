CFLAGS=-DDUMP=1 -Wall -O2 -events -g
CLIBS=-lm

Tin=473.15
Tam=300
Tfil_out=330
maxlevel=9
iter_fp=0
ratio_Rbmin=0.166666
ratio_Rbmax=1.
Nb=0
TOLERANCE_P=1e-7
TOLERANCE_V=1e-6
TOLERANCE_T=1e-7
parallel=4

args=$(Tin) $(Tam) $(Tfil_out) $(maxlevel) $(iter_fp) $(ratio_Rbmin) $(ratio_Rbmax) \
$(Nb) $(TOLERANCE_P) $(TOLERANCE_V) $(TOLERANCE_T)

kompozit:
	rm -rf .qcc* a.out dump* out log
#	qcc $(CFLAGS) poly_heat.c $(CINCL) $(CLIBS) &&	./a.out	$(args)	> out 2> log
	CC99='mpicc -std=c99' qcc -DDUMP=1 -Wall -O2 -events -D_MPI=1  kompozit.c -lm
	mpirun -np $(parallel) ./a.out $(args) >out 2>log
kompozit_test:
	#rm -rf .qcc* a.out dump* out log
	CC99='mpicc -std=c99' qcc -DDUMP=1 -Wall -O2 -events -g -D_MPI=1  kompozit.c -lm
	mpirun -np 2 ./a.out $(args) >out_test 2>log_test
	#mpirun -np 2 xterm -hold -e gdb -ex run --args ./a.out $(args)
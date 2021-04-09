CFLAGS=-DDUMP=1 -Wall -O2 -events -g
#OPENGLIBS=-L$(BASILISK)/gl -L/opt/local/lib/ -lfb_osmesa -lOSMesa -lglutils -framework OpenGL -framework GLUT
CLIBS=-lm $(OPENGLIBS)
CINC=-I/opt/local/include
args=
parallel=8


fixed_eb:
	rm -rf fixed .qcc*
	CC99='mpicc -std=c99' qcc $(CINC) $(CFLAGS) -D_MPI=1  fixed_eb.c $(CLIBS) -o fixed_eb.out
	mpirun -np $(parallel) ./fixed_eb.out $(args) >out_eb 2>log_eb

fixed_bp:
	rm -rf fixed .qcc*
	CC99='mpicc -std=c99' qcc $(CINC) $(CFLAGS) -D_MPI=1  fixed_bp.c $(CLIBS) -o fixed_bp.out
	mpirun -np $(parallel) ./fixed_bp.out $(args) >out_bp 2>log_bp

fixed_bp_test:
	rm -rf fixed .qcc*
	CC99='mpicc -std=c99' qcc $(CINC) $(CFLAGS)  -g -D_MPI=1  fixed_bp.c $(CLIBS) -o fixed_bp_test.out
	mpirun -np 2 ./fixed_bp_test.out $(args)
	#>out_bp 2>log_bp

fixed_bp_test_single:
	rm -rf fixed .qcc*
	qcc $(CINC) $(CFLAGS)  -g fixed_bp.c $(CLIBS) -o fixed_bp_single.out
	gdb ./fixed_bp_single.out $(args)
	#>out_bp 2>log_bp
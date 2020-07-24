CFLAG=-Wno-logical-op-parentheses -Wall -O2 -events -g
#-catch -g
all:
	echo "choose decaying_vortex_EB or decaying_vortex_Popinet or decaying_vortex_BP"
duct:
	rm -rf .qcc* log out res/duct_* snapshot-*
	qcc  $(CFLAG) run2.c -lm $(OPENGLIBS)
	./a.out > out 2> log
#	CC99='mpicc -std=c99' $(QCC_COMPILER)  -Wunused-result -O2 -D_MPI=1  cylinders_Popinet_trick.c -lm -L$(BASILISK)/gl $(OPENGLIB)
#	mpirun -np 4 ./a.out > out 2> log
duct_bp:
	rm -rf .qcc* log_bp out_bp res/duct_bp_* snapshot-*
	qcc  $(CFLAG) duct_bp.c -lm $(OPENGLIBS)
	./a.out > out_bp_dt=1e-3 2> log_bp_dt=1e-3
#	CC99='mpicc -std=c99' $(QCC_COMPILER)  -Wunused-result -O2 -D_MPI=1  cylinders_Popinet_trick.c -lm -L$(BASILISK)/gl $(OPENGLIB)
#	mpirun -np 4 ./a.out > out 2> log
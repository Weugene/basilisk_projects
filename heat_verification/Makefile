CFLAGS=-DDUMP=1 -Wall -O2 -events -g
CLIBS=-lm

maxlevel=7
iter_fp=0
TOLERANCE_P=1e-6
TOLERANCE_V=1e-7
TOLERANCE_T=1e-8
Htr=355000
Arrhenius_const=80600
#Ea_by_R=7697
Ea_by_R=6000
heat_verification_uniform_in_space:
	rm -rf .qcc* a.out dump* out log
#	qcc $(CFLAGS) heat_verification_uniform_in_space.c $(CINCL) $(CLIBS)
#	./a.out	 $(Tcyl) $(maxlevel) $(iter_fp) $(ratio_Rbmin) $(ratio_Rbmax) $(ratio_dist_x) $(ratio_dist_y) $(ratio_front_x) $(Nb) $(Ncx) $(Ncy) $(TOLERANCE) $(Htr) $(Arrhenius_const) $(Ea_by_R)
#	> out 2> log
	CC99='mpicc -std=c99' qcc -DDUMP=1 -Wall -O2 -events -D_MPI=1  heat_verification_uniform_in_space.c -lm
	mpirun -np 2 ./a.out $(maxlevel) $(iter_fp) $(TOLERANCE_P) $(TOLERANCE_V) $(TOLERANCE_T) $(Htr) $(Arrhenius_const) $(Ea_by_R) >out 2>log
heat_verification_nonuniform_in_space:
	rm -rf .qcc* a.out dump-* out log
	CC99='mpicc -std=c99' qcc -DDUMP=1 -Wall -O2 -events -D_MPI=1  heat_verification_nonuniform_in_space.c -lm
	mpirun -np 2 ./a.out $(maxlevel) $(iter_fp) $(TOLERANCE_P) $(TOLERANCE_V) $(TOLERANCE_T) $(Htr) $(Arrhenius_const) $(Ea_by_R) >out 2>log

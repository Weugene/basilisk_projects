CFLAGS=-DDUMP=1 -Wall -O2 -events -g
CLIBS=-lm

Tcyl=380
maxlevel=8
iter_fp=0
ratio_Rbmin=0.166666
ratio_Rbmax=1.
ratio_dist_x=2
ratio_dist_y=2
ratio_front_x=-3.5#-3.5
Nb=8#8
Ncx=5#5
Ncy=7#7
TOLERANCE_P=1e-6
TOLERANCE_V=1e-7
TOLERANCE_T=1e-6
Htr=355000
Arrhenius_const=80600
Ea_by_R=7697

heat_polymerization:
	rm -rf .qcc* a.out dump* out log
#	qcc $(CFLAGS) poly_heat.c $(CINCL) $(CLIBS)
#	./a.out	 $(Tcyl) $(maxlevel) $(iter_fp) $(ratio_Rbmin) $(ratio_Rbmax) $(ratio_dist_x) $(ratio_dist_y) $(ratio_front_x) $(Nb) $(Ncx) $(Ncy) $(TOLERANCE) $(Htr) $(Arrhenius_const) $(Ea_by_R)
#	> out 2> log
	CC99='mpicc -std=c99' qcc -DDUMP=1 -Wall -O2 -events -D_MPI=1  poly_heat.c -lm
	mpirun -np 6 ./a.out $(Tcyl) $(maxlevel) $(iter_fp) $(ratio_Rbmin) $(ratio_Rbmax) $(ratio_dist_x) $(ratio_dist_y) $(ratio_front_x) $(Nb) $(Ncx) $(Ncy) $(TOLERANCE_P) $(TOLERANCE_V) $(TOLERANCE_T) $(Htr) $(Arrhenius_const) $(Ea_by_R) >out 2>log

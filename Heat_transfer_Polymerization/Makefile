CFLAGS=-DDUMP=1 -Wall -O2 -events -g
CLIBS=-lm

Tcyl=380
Tin=300
maxlevel=9
iter_fp=0
ratio_Rbmin=0.166666
ratio_Rbmax=1.
ratio_dist_x=2
ratio_dist_y=2
ratio_front_x=-6
cyl_x= -5
Nb=0
Ncx=5
Ncy=8
TOLERANCE_P=1e-7
TOLERANCE_V=1e-6
TOLERANCE_T=1e-7
Htr=355000
Arrhenius_const=80600
Ea_by_R=7697 #already divided by 8.31
dev_r=0.
dev_x=0.
dev_y=0.
shift_x=0
shift_y=$(shell echo 0.5\*$(ratio_dist_y) | bc)
non_saturated=1

parallel=4

args=$(Tcyl) $(Tin) $(maxlevel) $(iter_fp) $(ratio_Rbmin) $(ratio_Rbmax) $(ratio_dist_x) $(ratio_dist_y) $(ratio_front_x) \
$(cyl_x) $(Nb) $(Ncx) $(Ncy) $(TOLERANCE_P) $(TOLERANCE_V) $(TOLERANCE_T) \
$(Htr) $(Arrhenius_const) $(Ea_by_R) $(dev_r) $(dev_x) $(dev_y) $(shift_x) $(shift_y) $(non_saturated)
heat_polymerization:
	rm -rf .qcc* a.out dump* out log
#	qcc $(CFLAGS) poly_heat.c $(CINCL) $(CLIBS) &&	./a.out	$(args)	> out 2> log
	CC99='mpicc -std=c99' qcc -DDUMP=1 -Wall -O2 -events -D_MPI=1  poly_heat.c -lm
	mpirun -np $(parallel) ./a.out $(args) >out 2>log
heat_polymerization_front:
	rm -rf .qcc* a.out dump* out log
	CC99='mpicc -std=c99' qcc -DDUMP=1 -Wall -O2 -events -D_MPI=1  poly_heat_front.c -lm
	mpirun -np $(parallel) ./a.out $(args) >out 2>log
heat_polymerization_front_test:
	#rm -rf .qcc* a.out dump* out log
	CC99='mpicc -std=c99' qcc -DDUMP=1 -Wall -O2 -events -g -D_MPI=1  poly_heat_front.c -lm
	mpirun -np 2 ./a.out $(args)
	#mpirun -np 2 xterm -hold -e gdb -ex run --args ./a.out $(args)
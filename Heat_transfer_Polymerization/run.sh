#!/bin/bash
set -x
set -o nounset
Tcyl=280
maxlevel=8
iter_fp=0
ratio_Rbmin=0.4
ratio_Rbmax=0.75
ratio_dist_x=2
ratio_dist_y=2
ratio_front_x=-3.5
Nb=8
Ncx=5
Ncy=7
TOLERANCE=1e-7
Htr=380.2e+3
Arrhenius_const=15063
Ea_by_R=1202.79
rm -rf .qcc* a.out dump* out log
#	qcc $(CFLAGS) poly_heat.c $(CINCL) $(CLIBS)
#	./a.out	 > out 2> log &
CC99='mpicc -std=c99' qcc -DDUMP=1 -Wall -O2 -events -D_MPI=1  poly_heat.c -lm
mpirun -np 2 ./a.out $Tcyl $maxlevel $iter_fp $ratio_Rbmin $ratio_Rbmax $ratio_dist_x $ratio_dist_y $ratio_front_x $Nb $Ncx $Ncy $TOLERANCE $Htr $Arrhenius_const $Ea_by_R >out 2>log

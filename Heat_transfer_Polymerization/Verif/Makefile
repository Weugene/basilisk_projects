CFLAGS=-DDUMP=1 -O2 #-events
CLIBS=-lm

Tcyl=400
Tin=300
maxlevel=8
iter_fp=0
ratio_Rbmin=0.166666
ratio_Rbmax=0.6
ratio_dist_x=2
ratio_dist_y=2
ratio_front_x=100
cyl_x=-0.5
Nb=5
Ncx=2
Ncy=3
TOLERANCE_P=1e-6
TOLERANCE_V=1e-6
TOLERANCE_T=1e-8
Htr=355000 #355000
Arrhenius_const=80600 #80600
Ea_by_R=2400 #already divided by 8.31  =  7697=64000/8.314
dev_r=0.
dev_x=0.
dev_y=0.
shift_x=0
shift_y=$(shell echo 0.5\*$(ratio_dist_y) | bc)
non_saturated=0

parallel=4
ind:=verif_wall_general1

dumpname=dump-1.05384e-05
args=$(Tcyl) $(Tin) $(maxlevel) $(iter_fp) $(ratio_Rbmin) $(ratio_Rbmax) $(ratio_dist_x) $(ratio_dist_y) $(ratio_front_x) \
$(cyl_x) $(Nb) $(Ncx) $(Ncy) $(TOLERANCE_P) $(TOLERANCE_V) $(TOLERANCE_T) \
$(Htr) $(Arrhenius_const) $(Ea_by_R) $(dev_r) $(dev_x) $(dev_y) $(shift_x) $(shift_y) $(non_saturated) $(ind)

args_short=$(Tcyl) $(Tin) $(maxlevel) $(iter_fp) $(TOLERANCE_P) $(TOLERANCE_V) $(TOLERANCE_T) \
$(Htr) $(Arrhenius_const) $(Ea_by_R) $(ind)

#run:=mpirun -np $(parallel)  ./a.out $(args) > out$(ind) 2> log$(ind)
run_short:=mpirun -np $(parallel)  ./a.out $(args_short) > out$(ind) 2> log$(ind)
compile:=CC99='mpicc -std=c99' qcc -O2 -D_MPI=1
remove:=rm -rf .qcc*

all: verify

verify:
	$(remove)
	$(compile)  verification_poly.c -lm
	$(run)
verify_wall:
	$(remove)
	$(compile)  verification_poly_wall.c -lm
	$(run)
verify_cylinder:
	$(remove)
	$(compile)  verification_poly_cylinder.c -lm
	$(run_short)
verify_couette:
	$(remove)
	$(compile)  verification_poly_couette.c -lm
	$(run_short)

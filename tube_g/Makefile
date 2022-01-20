CFLAG= -O2 -events
ACTIVE_ENV=source ~/.bash_profile
maxlevel=9
Re=0.005
Ca=0.04
Bo=0.0121
RhoR=1000
MuR=50
deq_st=0.3
lDomain=20
adapt_method=0
iter_fp=0
dt_vtk=0.1
snapshot_i=1000
parallel=8
ind=Glycerol
args:=$(maxlevel) $(Re) $(Ca) $(Bo) $(RhoR) $(MuR) $(deq_st) $(lDomain) $(adapt_method) $(iter_fp) $(dt_vtk) $(snapshot_i) $(ind)

run:=mpirun -np $(parallel)  ./a.out $(args) > out$(ind) 2>> log$(ind)
compile:=CC99='mpicc -std=c99' qcc -O2 -D_MPI=1 -Wall -events
remove:=rm -rf .qcc*

tube_bp_grav:
	$(remove)
	$(compile) tube_bp_g.c -lm
	$(run)

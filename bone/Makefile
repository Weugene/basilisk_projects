CFLAGS=-DDUMP=1 -Wall -O2 -events -g
CLIBS=-lm $(OPENGLIBS)
CINCL=-I$(BASILISK)
CORES=2
CFILE=bone
bone3D:
	rm -r .qcc* bone.pvd res/bone* log out || continue
	qcc   $(CFLAGS) bone.c $(CINCL) $(CLIBS)
	./a.out
#> out 2> log
#	CC99='mpicc -std=c99' qcc $(CFLAGS)  -D_MPI=1  $(CFILE).c $(CINCL) $(CLIBS)
#	mpirun -np $(CORES) ./a.out > out 2> log
bone3D_simple:
	rm -r .qcc* simple_bone.pvd res/simple_bone* log out || continue
	qcc   $(CFLAGS) bone_simple.c $(CINCL) $(CLIBS)
	./a.out

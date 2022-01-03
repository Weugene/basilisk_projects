CFILE=karman_GI_bubble
Re=100
Ca=0.1
Rrad=1
Rrho=1
Rmu=1
maxlevel=9
TOLERANCE=1e-4
ind=""
karman_GI_bubble:
	rm -rf .qcc* log out rk_* snapshot-* nameYouWant*
	qcc -O2 -events $(CFILE).c -lm $(OPENGLIBS)
	./a.out $(Re) $(Ca) $(Rrad) $(Rrho) $(Rmu) $(maxlevel) $(TOLERANCE)> out_$(ind) 2> log_$(ind)
#	CC99='mpicc -std=c99' qcc -events -Wunused-result -O2 -D_MPI=1  $(CFILE).c -lm $(OPENGLIBS)
#	mpirun -np 4 ./a.out > out 2> log

#	qcc -O2 -events karman_GI_bubble.c -lm $OPENGLIBS
#	CC99='mpicc -std=c99' qcc -events -Wunused-result -O2 -D_MPI=1  karman_GI_bubble.c -lm $OPENGLIBS

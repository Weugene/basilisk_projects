CFILE=moving_ring
Re=100
Ca=0.01
rad=0.0625
Rratio=2.0
Rrho=1
Rmu=1
maxlevel=9
lubrication=1.05
parallel=8
target=$(MAKECMDGOALS)
ind:=$(target)_jmax=$(maxlevel)_Re=$(Re)_Ca=$(Ca)_Rratio=$(Rratio)_lubr=$(lubrication)_no_patrcurv

args=$(Re) $(Ca) $(rad) $(Rratio) $(Rrho) $(Rmu) $(maxlevel) $(ind) $(lubrication)
run:=mpirun -np $(parallel)  ./a.out $(args) > out$(ind) 2>> log$(ind)
compile:=CC99='mpicc -std=c99' qcc -O2 -D_MPI=1 -Wall -events
remove:=rm -rf .qcc*
bpmufsticky:
	$(remove)
	$(compile) -DCORRECT_UF_FLUXES=1 -DSTICKY_SOLID=1 $(CFILE)_bpm.c -lm
	$(run)
bpmufnosticky:
	$(remove)
	$(compile) -DCORRECT_UF_FLUXES=1 -DNO_STICKY_SOLID=1 $(CFILE)_bpm.c -lm
	$(run)
bpmuflinear:
	$(remove)
	$(compile) -DCORRECT_UF_FLUXES=1 $(CFILE)_bpm.c -lm
	$(run)

bpmnoufcorr:
	$(remove)
	$(compile) -DCORRECT_UF_FLUXES=0 $(CFILE)_bpm.c -lm
	$(run)

popinetcooruf:
	$(remove)
	$(compile) -DCORRECT_UF_FLUXES=1 $(CFILE)_popinet.c -lm
	$(run)
popinetnocooruf:
	$(remove)
	$(compile) -DCORRECT_UF_FLUXES=0 $(CFILE)_popinet.c -lm
	$(run)

solidinbubble3D:
	$(remove)
	$(compile) -grid=octree $(CFILE)_bpm.c -lm
	$(run)
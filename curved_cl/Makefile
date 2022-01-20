CFILE=contact
Re=100
Ca=0.01
rad=0.2
Rratio=2.0
Rrho=1
Rmu=1
maxlevel=7
lubrication=1.05
parallel=8
target=$(MAKECMDGOALS)
ind:=$(target)_jmax=$(maxlevel)_Re=$(Re)_Ca=$(Ca)_Rratio=$(Rratio)_lubr=$(lubrication)_no_patrcurv_experiment_130

args=$(Re) $(Ca) $(rad) $(Rrho) $(Rmu) $(maxlevel) $(ind) $(lubrication)
run:=mpirun -np $(parallel)  ./a.out $(args) > out$(ind) 2>> log$(ind)
compile:=CC99='mpicc -std=c99' qcc -O2 -D_MPI=1 -Wall -events
remove:=rm -rf .qcc*
contactbpmufsticky_not_edge:
	$(remove)
	$(compile) -DCORRECT_UF_FLUXES=1 -DSTICKY_SOLID=1 $(CFILE).c -lm
	$(run)
contacbpmufnosticky:
	$(remove)
	$(compile) -DCORRECT_UF_FLUXES=1 -DNO_STICKY_SOLID=1 $(CFILE).c -lm
	$(run)
contacbpmuflinear:
	$(remove)
	$(compile) -DCORRECT_UF_FLUXES=1 $(CFILE).c -lm
	$(run)

contacbpmnoufcorr:
	$(remove)
	$(compile) -DCORRECT_UF_FLUXES=0 $(CFILE).c -lm
	$(run)

contacpopinetcooruf:
	$(remove)
	$(compile) -DCORRECT_UF_FLUXES=1 $(CFILE)_popinet.c -lm
	$(run)
contacpopinetnocooruf:
	$(remove)
	$(compile) -DCORRECT_UF_FLUXES=0 $(CFILE)_popinet.c -lm
	$(run)
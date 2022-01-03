CFILE=cyl
Re=100
Ca=0.1
rad=0.0625
Rrho=0.1
Rmu=0.1
maxlevel=9
parallel=8
target=$(MAKECMDGOALS)
ind:=$(target)_jmax=$(maxlevel)_Re=$(Re)_Ca=$(Ca)_no_patrcurv_film
args:=$(Re) $(Ca) $(rad) $(Rrho) $(Rmu) $(maxlevel) $(ind)
run:=mpirun -np $(parallel)  ./a.out $(args) > out$(ind) 2>> log$(ind)
compile:=CC99='mpicc -std=c99' qcc -O2 -D_MPI=1 -Wall -events
remove:=rm -rf .qcc*

cylonwaterPHOBMOVEnoufcorr:
	$(remove)
	$(compile) -DHYDROPHOBIC=1 -DCORRECT_UF_FLUXES=0 -DMOVE=1 $(CFILE).c -lm
	$(run)
cylonwaterPHOBSTEADYnoufcorr:
	$(remove)
	$(compile) -DHYDROPHOBIC=1 -DCORRECT_UF_FLUXES=0 -DMOVE=0 $(CFILE).c -lm
	$(run)
cylonwaterPHOBMOVEuflinear:
	$(remove)
	$(compile) -DHYDROPHOBIC=1 -DCORRECT_UF_FLUXES=1 -DMOVE=1 $(CFILE).c -lm
	$(run)
cylonwaterPHOBSTEADYuflinear:
	$(remove)
	$(compile) -DHYDROPHOBIC=1 -DCORRECT_UF_FLUXES=1 -DMOVE=0 $(CFILE).c -lm
	$(run)
cylonwaterPHOBMOVEufsticky:
	$(remove)
	$(compile) -DHYDROPHOBIC=1 -DCORRECT_UF_FLUXES=1 -DSTICKY_SOLID=1 -DMOVE=1 $(CFILE).c -lm
	$(run)
cylonwaterPHOBSTEADYufsticky:
	$(remove)
	$(compile) -DHYDROPHOBIC=1 -DCORRECT_UF_FLUXES=1 -DSTICKY_SOLID=1 -DMOVE=0 $(CFILE).c -lm
	$(run)

cylonwater3DPHOBMOVEnoufcorr:
	$(eval maxlevel=8)
	$(remove)
	$(compile) -DHYDROPHOBIC=1 -DMOVE=1 -DCORRECT_UF_FLUXES=0 -grid=octree $(CFILE).c -lm
	$(run)
cylonwater3DPHOBSTEADYnoufcorr:
	$(eval maxlevel=8)
	$(remove)
	$(compile) -DHYDROPHOBIC=1 -DMOVE=0 -DCORRECT_UF_FLUXES=0 -grid=octree $(CFILE).c -lm
	$(run)
cylonwater3DPHOBMOVEuflinear:
	$(eval maxlevel=8)
	$(remove)
	$(compile) -DHYDROPHOBIC=1 -DMOVE=1 -DCORRECT_UF_FLUXES=1 -grid=octree $(CFILE).c -lm
	$(run)
cylonwater3DPHOBSTEADYuflinear:
	$(eval maxlevel=8)
	$(remove)
	$(compile) -DHYDROPHOBIC=1 -DMOVE=0 -DCORRECT_UF_FLUXES=1 -grid=octree $(CFILE).c -lm
	$(run)
cylonwater3DPHOBMOVEufsticky:
	$(eval maxlevel=8)
	$(remove)
	$(compile) -DHYDROPHOBIC=1 -DMOVE=1 -DCORRECT_UF_FLUXES=1 -DSTICKY_SOLID=1 -grid=octree $(CFILE).c -lm
	$(run)
cylonwater3DPHOBSTEADYufsticky:
	$(eval maxlevel=8)
	$(remove)
	$(compile) -DHYDROPHOBIC=1 -DMOVE=0 -DCORRECT_UF_FLUXES=1 -DSTICKY_SOLID=1 -grid=octree $(CFILE).c -lm
	$(run)

cylonwaterPHILMOVE:
	$(remove)
	$(compile) -DHYDROPHOBIC=0 -DMOVE=1 $(CFILE).c -lm
	$(run)
cylonwaterPHILSTEADY:
	$(remove)
	$(compile) -DHYDROPHOBIC=0 -DMOVE=0 $(CFILE).c -lm
	$(run)
cylonwater3DPHILMOVE:
	$(eval maxlevel=8)
	$(remove)
	$(compile) -DHYDROPHOBIC=0 -DMOVE=1 -grid=octree $(CFILE).c -lm
	$(run)
cylonwater3DPHILSTEADY:
	$(eval maxlevel=8)
	$(remove)
	$(compile) -DHYDROPHOBIC=0 -DMOVE=0 -grid=octree $(CFILE).c -lm
	$(run)

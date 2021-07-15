CFILE=cyl
Re=100
Ca=0.1
rad=0.0625
Rrho=0.1
Rmu=0.1
maxlevel=10
parallel=8
ind:=$(maxlevel)_mod_kappa4
cylonwaterPHILMOVE:
	rm -rf .qcc*
	CC99='mpicc -std=c99' qcc -O2 -D_MPI=1 -DHYDROPHOBIC=0 -DMOVE=1  -Wall -events $(CFILE).c -lm
	mpirun -np $(parallel)  ./a.out $(Re) $(Ca) $(rad) $(Rrho) $(Rmu) $(maxlevel) "vtk_2DPHILMOVE_$(ind)"> out2DPHILMOVE_$(ind) 2> log2DPHILMOVE_$(ind)
cylonwaterPHOBMOVE:
	rm -rf .qcc*
	CC99='mpicc -std=c99' qcc -O2 -D_MPI=1 -DHYDROPHOBIC=1 -DMOVE=1 -Wall -events $(CFILE).c -lm
	mpirun -np $(parallel)  ./a.out $(Re) $(Ca) $(rad) $(Rrho) $(Rmu) $(maxlevel) "vtk_2DPHOBMOVE_$(ind)"> out2DPHOBMOVE_$(ind) 2> log2DPHOBMOVE_$(ind)
cylonwaterPHILSTEADY:
	rm -rf .qcc*
	CC99='mpicc -std=c99' qcc -O2 -D_MPI=1 -DHYDROPHOBIC=0 -DMOVE=0 -Wall -events $(CFILE).c -lm
	mpirun -np $(parallel)  ./a.out $(Re) $(Ca) $(rad) $(Rrho) $(Rmu) $(maxlevel) "vtk_2DPHILSTEADY_$(ind)"> out2DPHILSTEADY_$(ind) 2> log2DPHILSTEADY_$(ind)
cylonwaterPHOBSTEADY:
	rm -rf .qcc*
	CC99='mpicc -std=c99' qcc -O2 -D_MPI=1 -DHYDROPHOBIC=1 -DMOVE=0 -Wall -events $(CFILE).c -lm
	mpirun -np $(parallel)  ./a.out $(Re) $(Ca) $(rad) $(Rrho) $(Rmu) $(maxlevel) "vtk_2DPHOBSTEADY_$(ind)"> out2DPHOBSTEADY_$(ind) 2> log2DPHOBSTEADY_$(ind)

cylonwater3DPHILMOVE:
	rm -rf .qcc*
	CC99='mpicc -std=c99' qcc -O2 -D_MPI=1 -DHYDROPHOBIC=0 -DMOVE=1 -Wall -events -grid=octree $(CFILE).c -lm
	mpirun -np  $(parallel) ./a.out $(Re) $(Ca) $(rad) $(Rrho) $(Rmu) $(maxlevel) "vtk_3DPHILMOVE_$(ind)" > out3DPHILMOVE_$(ind) 2> log3DPHILMOVE_$(ind)
cylonwater3DPHOBMOVE:
	rm -rf .qcc*
	CC99='mpicc -std=c99' qcc -O2 -D_MPI=1 -DHYDROPHOBIC=1 -DMOVE=1 -Wall -events -grid=octree $(CFILE).c -lm
	mpirun -np $(parallel)  ./a.out $(Re) $(Ca) $(rad) $(Rrho) $(Rmu) $(maxlevel) "vtk_3DPHOBMOVE_$(ind)" > out3DPHOBMOVE_$(ind) 2> log3DPHOBMOVE_$(ind)
cylonwater3DPHILSTEADY:
	rm -rf .qcc*
	CC99='mpicc -std=c99' qcc -O2 -D_MPI=1 -DHYDROPHOBIC=0 -DMOVE=0  -Wall -events -grid=octree $(CFILE).c -lm
	mpirun -np  $(parallel) ./a.out $(Re) $(Ca) $(rad) $(Rrho) $(Rmu) $(maxlevel) "vtk_3DPHILSTEADY_$(ind)" > out3DPHILSTEADY_$(ind) 2> log3DPHILSTEADY_$(ind)
cylonwater3DPHOBSTEADY:
	rm -rf .qcc*
	CC99='mpicc -std=c99' qcc -O2 -D_MPI=1 -DHYDROPHOBIC=1 -DMOVE=0 -Wall -events -grid=octree $(CFILE).c -lm
	mpirun -np $(parallel)  ./a.out $(Re) $(Ca) $(rad) $(Rrho) $(Rmu) $(maxlevel) "vtk_3DPHOBSTEADY_$(ind)" > out3DPHOBSTEADY_$(ind) 2> log3DPHOBSTEADY_$(ind)

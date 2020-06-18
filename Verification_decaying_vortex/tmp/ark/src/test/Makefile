all: check 3D-tests mpi-tests load-balancing curvature-tests axi-tests \
	cadna-tests embed-tests multilayer-tests gotm-tests

include $(BASILISK)/Makefile.defs

# the default CFLAGS are set in $(BASILISK)/config
CFLAGS += -O2 -DMTRACE=3

# list only non-default tests
check: halo.vtst interpolate.vtst events.vtst faces.vtst coarsen.vtst \
	lid-restore.tst poiseuille.vtst tag.ctst cvmix.vtst

# these tests have special dependencies/compilation requirements

basilisk.tst: basilisk.gnu
basilisk.tst: CFLAGS += -DDEBUG=1

dry-explicit.c: dry.c
	ln -s dry.c dry-explicit.c 
dry-explicit.tst: dry.s
dry-explicit.tst: CFLAGS += -DEXPLICIT=1

dry.tst: dry-explicit.tst

drybump-explicit.c: drybump.c
	ln -s drybump.c drybump-explicit.c 
drybump-explicit.tst: drybump.s
drybump-explicit.tst: CFLAGS += -DEXPLICIT=1

drybump.tst: drybump-explicit.tst

drybump2D-implicit.c: drybump2D.c
	ln -s drybump2D.c drybump2D-implicit.c 
drybump2D-implicit.tst: drybump2D.s
drybump2D-implicit.tst: CFLAGS += -DIMPLICIT=1

drybump2D.tst: drybump2D-implicit.tst

lake-tr.tst: lake-tr-explicit.tst

lake-tr-explicit.c: lake-tr.c
	ln -s lake-tr.c lake-tr-explicit.c
lake-tr-explicit.tst: lake-tr.s
lake-tr-explicit.tst: CFLAGS += -DEXPLICIT=1

multiriverinflow.tst: multiriverinflow.ctst
multiriverinflow.tst: CFLAGS += -fopenmp
multiriverinflow.ctst: CFLAGS += -fopenmp

parabola-explicit.c: parabola.c
	ln -s parabola.c parabola-explicit.c 
parabola-explicit.tst: parabola.s
parabola-explicit.tst: CFLAGS += -DEXPLICIT=1

parabola.tst: parabola-explicit.tst parabola-ml.tst

seawallsv.c: seawall.c
	ln -s seawall.c seawallsv.c 
seawallsv.tst: seawall.s
seawallsv.tst: CFLAGS += -DSAINT_VENANT=1

seawall.tst: seawallsv.tst seawall-ml.tst

conicalsv.c: conical.c
	ln -s conical.c conicalsv.c 
conicalsv.tst: conical.s
conicalsv.tst: CFLAGS += -DSAINT_VENANT=1

conical-implicit.c: conical.c
	ln -s conical.c conical-implicit.c 
conical-implicit.tst: conical.s
conical-implicit.tst: CFLAGS += -DSAINT_VENANT=1 -DIMPLICIT=1

conical.tst: conicalsv.tst conical-implicit.tst conical-ml.tst

curvature.tst: curvature.3D.tst

explicit.c: implicit.c
	ln -s implicit.c explicit.c 
explicit.tst: implicit.s
explicit.tst: CFLAGS += -DEXPLICIT=1

implicit.tst: explicit.tst

bore.tst: bore1.tst

explosion3D.tst: CFLAGS=-grid=multigrid3D

explosion3D.c: explosion.c
	ln -s explosion.c explosion3D.c

explosion.tst: explosion.ctst explosion.3D.tst explosion3D.tst

laplacian.tst: laplacian.ctst

lidmac.c: lid.c
	ln -s lid.c lidmac.c
lidmac.tst: lid.s
lidmac.tst: CFLAGS += -DMAC=1

lid.tst: lidmac.tst

lid-restore.c: lid.c
	ln -s lid.c lid-restore.c
lid-restore.dump: lid.tst
	rm -f lid-restore.dump; ln -s lid/dump lid-restore.dump
lid-restore.tst: lid-restore.dump

oscillation-momentum.c: oscillation.c
	ln -s oscillation.c oscillation-momentum.c 
oscillation-momentum.s: oscillation.s
oscillation-momentum.tst: CFLAGS += -DMOMENTUM=1

oscillation.tst: oscillation-momentum.tst

ponds-implicit.c: ponds.c
	ln -s ponds.c ponds-implicit.c
ponds-implicit.s: ponds.s
ponds-implicit.tst: CFLAGS += -DIMPLICIT=1

ponds.tst: ponds-implicit.tst

poiseuille-fenep.c: poiseuille-oldroydb.c
	ln -s poiseuille-oldroydb.c poiseuille-fenep.c
poiseuille-fenep.s: poiseuille-oldroydb.s
poiseuille-fenep.tst: CFLAGS += -DFENE_P=1

poiseuille-oldroydb.tst: poiseuille-fenep.tst

poisson.tst: poisson.ctst

rising-axi.c: rising.c
	ln -s rising.c rising-axi.c
rising-axi.tst: rising.s
rising-axi.tst: CFLAGS += -DAXIS=1

rising2.c: rising.c
	ln -s rising.c rising2.c
rising2.tst: rising.s
rising2.tst: CFLAGS += -DCASE2=1

rising-reduced.c: rising.c
	ln -s rising.c rising-reduced.c
rising-reduced.tst: rising.s
rising-reduced.tst: CFLAGS += -DREDUCED=1

reversed.tst: reversed.ctst

rotate.tst: rotate.ctst

rt-reduced.c: rt.c
	ln -s rt.c rt-reduced.c
rt-reduced.tst: rt.s
rt-reduced.tst: CFLAGS += -DREDUCED=1

taylor-green-all-mach.c: taylor-green.c
	ln -s taylor-green.c taylor-green-all-mach.c 
taylor-green-all-mach.tst: taylor-green-all-mach.s
taylor-green-all-mach.tst: CFLAGS += -DALL_MACH=1

taylor-green.tst: taylor-green-all-mach.tst

view.tst: CC = mpicc -D_MPI=4
view.tst: view.3D.tst

view.3D.tst: CC = mpicc -D_MPI=4

coarsen.vtst: CFLAGS = -DMTRACE=3 -progress

# Axisymmetric tests

axi-tests: axiadvection.tst axi.tst poiseuille-axi.tst \
	rising-axi.tst

# 3D tests

3D-tests: circle.3D.tst curvature.3D.tst hf.3D.tst periodic.3D.tst \
	poisson.3D.tst refineu.3D.tst solenoidal.3D.tst

# Curvature tests

curvature-tests: hf.tst hf3D.tst hf.ctst hf.3D.tst hf1.tst hf-mask.tst \
	curvature.tst curvature.3D.tst

hf-mask.tst: CFLAGS = -DTRASH=1

# Multilayer tests

multilayer-tests: dispersion-gn.tst dispersion.tst \
	beach-ml.tst seawall-ml.tst conical-ml.tst \
	bar-ml.tst parabola-ml.tst lonlat-ml.tst \
	wind-driven-nh.tst wind-driven.tst large.tst \
	gaussian-stvt.tst gaussian-hydro.tst gaussian.tst \
	stokes-ns.tst stokes.tst ponds-ml.tst

dispersion-gn.c: dispersion.c
	ln -s dispersion.c dispersion-gn.c 
dispersion-gn.tst: dispersion-gn.s

dispersion.tst: CFLAGS += -DLAYERS=1

bar-ml.c: bar.c
	ln -s bar.c bar-ml.c 
bar-ml.tst: bar-ml.s
bar-ml.tst: CFLAGS += -DLAYERS=1

wind-driven-nh.c: wind-driven.c
	ln -s wind-driven.c wind-driven-nh.c 
wind-driven-nh.tst: wind-driven-nh.s
wind-driven-nh.tst: CFLAGS += -DLAYERS=1 -DNH=1

wind-driven-stvt.c: wind-driven.c
	ln -s wind-driven.c wind-driven-stvt.c 
wind-driven-stvt.tst: wind-driven-stvt.s

wind-driven.tst: CFLAGS += -DLAYERS=1

large.tst: CFLAGS += -DLAYERS=1
large.tst: large-ns.tst

gaussian-hydro.c: gaussian.c
	ln -s gaussian.c gaussian-hydro.c 
gaussian-hydro.tst: gaussian-hydro.s
gaussian-hydro.tst: CFLAGS += -DLAYERS=1 -DHYDRO=1

gaussian-stvt.c: gaussian.c
	ln -s gaussian.c gaussian-stvt.c 
gaussian-stvt.tst: gaussian-stvt.s

gaussian.tst: CFLAGS += -DLAYERS=1

parabola-ml.c: parabola.c
	ln -s parabola.c parabola-ml.c 
parabola-ml.tst: parabola-ml.s
parabola-ml.tst: CFLAGS += -DLAYERS=1

seawall-ml.c: seawall.c
	ln -s seawall.c seawall-ml.c 
seawall-ml.tst: seawall.s
seawall-ml.tst: CFLAGS += -DLAYERS=1

beach-ml.c: beach.c
	ln -s beach.c beach-ml.c 
beach-ml.tst: beach.s
beach-ml.tst: CFLAGS += -DLAYERS=1

conical-ml.c: conical.c
	ln -s conical.c conical-ml.c 
conical-ml.tst: conical.s
conical-ml.tst: CFLAGS += -DLAYERS=1

lonlat-ml.c: lonlat.c
	ln -s lonlat.c lonlat-ml.c 
lonlat-ml.tst: lonlat.s
lonlat-ml.tst: CFLAGS += -DLAYERS=1

ponds-ml.c: ponds.c
	ln -s ponds.c ponds-ml.c 
ponds-ml.tst: ponds.s
ponds-ml.tst: CFLAGS += -DLAYERS=1

stokes.tst: CFLAGS += -DLAYERS=1

## MPI tests

# bump2Dp.tst crashes with -DDEBUGCOND=1 but works with -DDEBUGCOND=0

mpi-tests: indexing.tst indexing.3D.tst \
	mpi-restriction.tst mpi-restriction.3D.tst \
	mpi-reduce.tst mpi-refine.tst mpi-refine1.tst mpi-refine.3D.tst \
	mpi-laplacian.tst mpi-laplacian.3D.tst \
	mpi-circle.tst mpi-circle1.tst mpi-flux.tst \
	mpi-interpu.tst mpi-coarsen.tst mpi-coarsen1.tst \
	hf1.tst pdump.tst restore.tst \
	pdump-multigrid.tst restore-multigrid.tst \
	restore-tree.tst \
	poiseuille-periodic.tst \
	gfsi.tst gfs.tst \
	load-balancing \
	mpi-grid.tst mpi-periodic-3D.tst \
	source.tst tag.tst tag1.tst view.tst view.3D.tst \
	boundary_vertex.tst boundary_vertex3D.tst \
	foreach_bnd1.tst

indexing.tst:		CC = mpicc -D_MPI=3
indexing.3D.tst:	CC = mpicc -D_MPI=3
mpi-restriction.tst:	CC = mpicc -D_MPI=3
mpi-restriction.3D.tst:	CC = mpicc -D_MPI=3
mpi-reduce.tst:		CC = mpicc -D_MPI=3
mpi-refine.tst:		CC = mpicc -D_MPI=4
mpi-refine1.tst:	CC = mpicc -D_MPI=11
mpi-refine.3D.tst:	CC = mpicc -D_MPI=4
mpi-laplacian.tst:	CC = mpicc -D_MPI=3
mpi-laplacian.3D.tst:	CC = mpicc -D_MPI=3
mpi-circle.tst:		CC = mpicc -D_MPI=5
foreach_bnd1.tst:	CC = mpicc -D_MPI=4

mpi-circle1.c: mpi-circle.c
	ln -s mpi-circle.c mpi-circle1.c

mpi-circle1.tst:  CC = mpicc -D_MPI=6
mpi-flux.tst:     CC = mpicc -D_MPI=6
mpi-interpu.tst:  CC = mpicc -D_MPI=5
mpi-coarsen.tst:  CC = mpicc -D_MPI=2
mpi-coarsen1.tst: CC = mpicc -D_MPI=5
bump2Dp.tst:      CC = mpicc -D_MPI=55
vortex.tst:	  CC = mpicc -D_MPI=7 -DJACOBI=1
axiadvection.tst: CC = mpicc -D_MPI=7 -DJACOBI=1
hf1.tst:	  CC = mpicc -D_MPI=7
poiseuille-periodic.tst: CC = mpicc -D_MPI=4
source.tst: CC = mpicc -D_MPI=7

tag.tst: CC = mpicc -D_MPI=7
tag1.tst: CC = mpicc -D_MPI=7

boundary_vertex.tst: CC = mpicc -D_MPI=7
boundary_vertex3D.tst: CC = mpicc -D_MPI=7

# parallel dump()/restore()

pdump.tst: CC = mpicc -D_MPI=7
restore.dump: pdump.tst
	rm -f restore.dump; ln -s pdump/restore.dump
restore.tst: restore.dump

pdump-multigrid.c: pdump.c
	ln -s pdump.c pdump-multigrid.c
pdump-multigrid.tst: pdump.c
pdump-multigrid.tst: CFLAGS = -grid=multigrid
pdump-multigrid.tst: CC = mpicc -D_MPI=16

restore-multigrid.dump: pdump-multigrid.tst
	rm -f restore-multigrid.dump; ln -s pdump-multigrid/restore.dump \
		restore-multigrid.dump
restore-multigrid.tst: restore-multigrid.dump
restore-multigrid.tst: CFLAGS = -grid=multigrid
restore-multigrid.tst: CC = mpicc -D_MPI=16

restore-tree.dump: pdump.tst
	rm -f restore-tree.dump; ln -s pdump/restore.dump restore-tree.dump
restore-tree.tst: restore-tree.dump
restore-tree.tst: CFLAGS = -DDEBUGCOND=false
restore-tree.tst: CC = mpicc -D_MPI=23

bump2Dp-restore.c: bump2Dp.c
	ln -s bump2Dp.c bump2Dp-restore.c
bump2Dp-restore.dump: bump2Dp.tst
	rm -f bump2Dp-restore.dump; ln -s bump2Dp/dump bump2Dp-restore.dump
bump2Dp-restore.tst: bump2Dp-restore.dump
bump2Dp-restore.tst: CC = mpicc -D_MPI=15

# parallel output_gfs()

gfs.tst: CC = mpicc -D_MPI=7
gfsi.gfs: gfs.tst
	rm -f gfsi.gfs; ln -s gfs/gfsi.gfs
gfsi.tst: gfsi.gfs

# load-balancing

load-balancing: balance5.tst balance6.tst balance7.tst \
		bump2Dp.tst bump2Dp-restore.tst vortex.tst axiadvection.tst

balance5.tst: CC = mpicc -D_MPI=9
balance6.c: balance5.c
	ln -s balance5.c balance6.c
balance6.tst: CC = mpicc -D_MPI=17
balance7.tst: CC = mpicc -D_MPI=17

# MPI-parallel multigrid

mpi-grid.tst: CC = mpicc -D_MPI=4
mpi-periodic-3D.c: periodic.c
	ln -s periodic.c mpi-periodic-3D.c
mpi-periodic-3D.s: periodic.s $(BASILISK)/grid/multigrid-mpi.h
mpi-periodic-3D.tst: CC = mpicc -D_MPI=8
mpi-periodic-3D.tst: CFLAGS += -grid=multigrid3D

# CADNA tests

cadna-tests: poiseuille.CADNA.tst sag.CADNA.tst poisson.CADNA.tst

# Periodic boundary tests

periodic1.tst: CFLAGS += -grid=bitree -DTRASH=1
periodic2.c: periodic1.c
	ln -s periodic1.c periodic2.c
periodic2.tst: CFLAGS += -DTRASH=1
periodic3.c: periodic1.c
	ln -s periodic1.c periodic3.c
periodic3.tst: CFLAGS += -grid=octree -DTRASH=1

# Embedded boundary tests

dirichlet.tst: CFLAGS += -DDIRICHLET=1
dirichlet.c: neumann.c
	ln -s neumann.c dirichlet.c

dirichlet3D.tst: CFLAGS += -DDIRICHLET=1
dirichlet3D.c: neumann3D.c
	ln -s neumann3D.c dirichlet3D.c

hydrostatic2.tst: CC = mpicc -D_MPI=8

embed-tests: neumann.tst dirichlet.tst poiseuille45.tst \
	couette.tst wannier.tst hydrostatic2.tst cylinders.tst \
	porous.tst porous1.tst uf.tst starting.tst \
	neumann3D.tst dirichlet3D.tst hydrostatic3.tst spheres.tst

# General Ocean Turbulence Model (GOTM) tests

gotm-tests: couette-gotm.tst entrainment.tst ows_papa.tst

couette-gotm.tst: CFLAGS += -DLAYERS=1

entrainment.tst: CFLAGS += -DLAYERS=1

ows_papa-kpp.c: ows_papa.c
	ln -s ows_papa.c ows_papa-kpp.c 
ows_papa-kpp.tst: ows_papa.s
ows_papa-kpp.tst: CFLAGS += -DLAYERS=1 -DKPP=1

ows_papa.tst: ows_papa-kpp.tst
ows_papa.tst: CFLAGS += -DLAYERS=1

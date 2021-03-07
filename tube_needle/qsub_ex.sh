#!/bin/bash
#PBS -N tube_needle
#PBS -q batch
#PBS -l nodes=1:ppn=25
#PBS -l pmem=10gb
#PBS -l walltime=40:00:00
###PBS -l feature=htc
case=22
maxlevel=9
Ca=0.023
Umean=1.580
Vd=0.2179e-9
adapt_method=0
iter_fp=0
lDomain=10
dt_vtk=1e-1
snapshot_i=500
echo "nodes:"
cat $PBS_NODEFILE
echo "openglibs=$OPENGLIBS"
cd $HOME/basilisk/work/tube_needle/
mpirun ./a.out $maxlevel $case $adapt_method $iter_fp $lDomain $dt_vtk $snapshot_i

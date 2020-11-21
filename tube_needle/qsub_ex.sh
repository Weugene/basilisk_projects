#!/bin/bash
#PBS -N tube
#PBS -q batch
#PBS -l nodes=2:ppn=2
#PBS -l pmem=1gb
#PBS -l walltime=00:10:00
###PBS -l feature=htc
case=22
maxlevel=8
Ca=0.023
Umean=1.580
Vd=0.2179e-9
adapt_method=1
iter_fp=0
lDomain=30
dt_vtk=1e-3
echo "nodes:"
cat $PBS_NODEFILE
echo "openglibs=$OPENGLIBS"
cd $HOME/basilisk/work/tube/
lscpu
mpirun ./a.out $maxlevel $Ca $Umean $Vd $adapt_method $iter_fp $lDomain $dt_vtk

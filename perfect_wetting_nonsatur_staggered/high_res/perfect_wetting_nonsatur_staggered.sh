#!/bin/bash

#PBS -l nodes=1:ppn=10

#PBS -l pmem=10000MB

#PBS -l walltime=40:00:00

#PBS -q batch

cd $HOME/basilisk/work/perfect_wetting_nonsatur_staggered/high_res
mpirun -np 10 ./a.out


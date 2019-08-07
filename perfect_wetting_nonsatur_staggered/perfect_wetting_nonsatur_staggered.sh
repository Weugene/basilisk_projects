#!/bin/bash

#PBS -l nodes=2:ppn=10

#PBS -l pmem=10000MB

#PBS -l walltime=40:00:00

#PBS -q batch

cd $HOME/basilisk/work/perfect_wetting_nonsatur_staggered
mpirun -np 20 ./a.out


#!/bin/bash

#PBS -l nodes=1:ppn=32

#PBS -l pmem=8000MB

#PBS -l walltime=23:00:00

#PBS -q batch
#PBS -e error.txt
#PBS -o out.txt
cd $HOME/basilisk/work/Kasimov_copy
mpirun -np 32 ./a.out


#!/bin/bash

#PBS -l nodes=1:ppn=10

#PBS -l pmem=8000MB

#PBS -l walltime=20:00:00

#PBS -q batch

cd $HOME/basilisk/work/Hele-Shaw_original/line/
mpirun -np 10 ./a.out


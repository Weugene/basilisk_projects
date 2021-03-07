#!/bin/bash

#PBS -l nodes=2:ppn=8

#PBS -l pmem=8000MB

#PBS -l walltime=48:00:00

#PBS -q batch

cd $HOME/basilisk/work/Hele-Shaw_original/line_high_res/
mpirun -np 16 ./a.out


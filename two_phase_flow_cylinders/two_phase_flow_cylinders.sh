#!/bin/bash

#PBS -l nodes=1:ppn=10

#PBS -l pmem=8000MB

#PBS -l walltime=45:00:00

#PBS -q batch

cd $HOME/basilisk/work/two_phase_flow_cylinders
mpirun -np 10 ./a.out


#!/bin/bash

#PBS -l nodes=1:ppn=10

#PBS -l pmem=8000MB

#PBS -l walltime=20:00:00

#PBS -q batch
#PBS -e error.txt
#PBS -o out.txt
cd ~/basilisk/work/Kasimov_Rayleigh_Plateau/less_sigma/
mpirun -np 10 ./a.out


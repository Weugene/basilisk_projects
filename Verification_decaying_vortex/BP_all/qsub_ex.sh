#!/bin/bash
#PBS -N vortex
#PBS -l nodes=1:ppn=15
#PBS -l walltime=30:00:00
cd $HOME/basilisk/work/Verification_decaying_vortex/BP
source $HOME/.bashrc
./runbash.sh
sleep 30h

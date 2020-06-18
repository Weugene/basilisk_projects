#!/bin/bash
#PBS -N vortex
#PBS -l nodes=1:ppn=16
#PBS -l walltime=30:00:00
cd $HOME/basilisk/work/Verification_decaying_vortex/BP8_modif_Ch_p_interp
source $HOME/.bashrc
./runbash.sh
sleep 30h

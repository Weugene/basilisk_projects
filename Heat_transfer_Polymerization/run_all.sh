#!/bin/bash
set -x
set -o nounset
echo "Usage: qsub or sbatch? shift_y non_saturated?"

mode=$1
if [[ "$mode" == "qsub" || "$mode" == "sbatch" ]]; then
  echo "You got the right input: $mode"
else
  echo "No way. Use qsub or sbatch"
  exit
fi

export PBS_O_WORKDIR=$PWD
echo "PBS_O_WORKDIR=${PBS_O_WORKDIR}"
export JLAB_SLURM_O_WORKDIR=$PWD
echo "JLAB_SLURM_O_WORKDIR=${JLAB_SLURM_O_WORKDIR}"

#./run_heat.sh $1 0 0
./run_heat.sh $1 0 1
#./run_heat.sh $1 1 0
./run_heat.sh $1 1 1



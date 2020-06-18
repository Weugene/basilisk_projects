#!/bin/bash
echo "Usage: No arguments"
dir_prefix="res_l"
#parameters for solver
solver_serial=./a.out
solver_parallel=./parallel.out
script=./run_one_experiment.sh
#levels=(4)
levels=()
levels_p=(11)
cells_per_zone=(1)
eps=(1e-2 5e-3 1e-3)
np_parallel=4
export PBS_O_WORKDIR=$PWD
#echo "OPENGLIBS=$OPENGLIBS"
echo "${PBS_O_WORKDIR}"
for e in "${eps[@]}"; do
for n in "${cells_per_zone[@]}"; do
for a in "${levels[@]}"; do
	tmp=${a}_${n}_${e}
	${script} ${PWD}/${solver_serial} serial 1 ${a} ${n} ${e} &
done
done
done

echo "parallel qsub running"
for e in "${eps[@]}"; do
for n in "${cells_per_zone[@]}"; do
for a in "${levels_p[@]}"; do
	tmp=${a}_${n}_${e}
	${script} ${PWD}/${solver_parallel} parallel ${np_parallel} ${a} ${n} ${e} &
done
done
done

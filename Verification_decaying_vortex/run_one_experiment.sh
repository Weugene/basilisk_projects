#!/bin/bash
source $HOME/.bashrc
echo "OPENGLIBS=$OPENGLIBS"
full_solver=$1
path=$(dirname ${full_solver})
cd $path
echo "my path=$PWD"
solver=./$(basename ${full_solver})
mode="$2"
np_parallel=$3
level=$4
cells_per_zone=$5
eps=$6
echo "full_solver=${full_solver} mode=${mode} np=$np_parallel level=$level cells_per_zone=$cells_per_zone eps=$eps"
tmp=${level}_${cells_per_zone}_${eps}
echo "${solver} arg=${tmp}"
dir="res_${tmp}"
mkdir "$dir" || continue
cp "${solver}" "$dir"
(
	cd "$dir" || exit
	if [[ "${mode}" == "serial" ]]
	then
		echo "serial mode: $solver "
		$solver $level ${cells_per_zone} $eps > out_${tmp} 2> log_${tmp}
	elif [[ "${mode}" == "parallel" ]]
	then
		echo "parallel mode: mpirun -np ${np_parallel} ${solver}"
		mpirun -np ${np_parallel} $solver ${level} ${cells_per_zone} ${eps} > out_${tmp} 2> log_${tmp}
	fi
)
